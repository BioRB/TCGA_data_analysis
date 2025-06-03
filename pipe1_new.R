#!/usr/bin/env Rscript

print("START")
# Imposta timeout globale a 120 secondi
options(timeout = 120)

# Configura timeout per connessioni HTTP
library(httr)
set_config(config(connecttimeout = 120, timeout = 120))

# Configurazione logging avanzata
log_dir <- "logs"
dir.create(log_dir, showWarnings = FALSE)

log_file <- file.path(log_dir, "analisi_tcga_log.txt")
warn_file <- file.path(log_dir, "analisi_tcga_warnings.txt")

# Inizializza i file
file.create(log_file, showWarnings = FALSE)
file.create(warn_file, showWarnings = FALSE)

# Chiudi eventuali sink aperti (versione robusta)
max_attempts <- 10
attempts <- 0
while (sink.number() > 0 && attempts < max_attempts) {
  sink()
  attempts <- attempts + 1
  Sys.sleep(0.1)
}

attempts <- 0
while (sink.number(type = "message") > 0 && attempts < max_attempts) {
  tryCatch({
    sink(type = "message")
  }, error = function(e) {
    message("Errore chiusura sink message: ", e$message)
  })
  attempts <- attempts + 1
  Sys.sleep(0.1)
}

# Apri i nuovi sink
warn_con <- file(warn_file, open = "a")  # Connessione per i warning
sink(log_file, append = TRUE, split = TRUE)
sink(warn_con, type = "message")

# Funzione per chiusura pulita alla fine
cleanup <- function() {
  try({
    # Chiudi tutti i sink standard
    while(sink.number() > 0) {
      sink()
    }
    
    # Chiudi tutti i sink di messaggi
    while(sink.number(type = "message") > 0) {
      sink(type = "message")
    }
    
    # Chiudi esplicitamente la connessione del file di warning
    if(exists("warn_con") && isOpen(warn_con)) {
      try(close(warn_con), silent = TRUE)
    }
  }, silent = TRUE)
}

# Registra la funzione di cleanup
reg.finalizer(environment(), function(e) cleanup(), onexit = TRUE)

# Carica librerie
library(readr)
library(stringr)
library(TCGAbiolinks)
library(dplyr)
library(maftools)
library(lubridate)
library(data.table)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(sjPlot)
library(GGally)
library(introdataviz)
library(ggExtra)
library(GDCRNATools)
library(tidyr)
library(showtext) # Aggiunto per gestione font

# Attiva il rendering dei font di sistema con showtext
showtext_auto()

message("Librerie caricate con successo.")
# Carica dati
oncogenic_genes_file <- "oncogenic_genes.txt" # Nome file di default
if (file.exists(oncogenic_genes_file)) {
  driver_genes <- readLines(oncogenic_genes_file)
  # Non c'è un dataframe qui, quindi non possiamo usare problems() direttamente.
  # Tuttavia, potresti aggiungere controlli sulla formattazione del file se necessario.
  message(paste("Caricati", length(driver_genes), "geni oncogenici dal file:", oncogenic_genes_file))
} else {
  stop(paste("Errore: il file dei geni oncogenici", oncogenic_genes_file, "non è stato trovato."))
}

# Variabile globale per controllo convergenza
ABILITA_CONVERGENZA <- FALSE

# Funzioni principali ------------------------------------------------------

run_regression <- function(data, gene, predictors) {
  formula <- as.formula(paste(gene, "~", paste(predictors, collapse = "+")))
  
  withCallingHandlers({
    tryCatch({
      model <- glm(formula, data = data, family = "binomial",
                   control = list(maxit = 25, epsilon = 1e-04))
      
      if (ABILITA_CONVERGENZA && !model$converged) {
        msg <- paste("Modello non convergente per il gene:", gene,
                     "con predittori:", paste(predictors, collapse = ","))
        warning(msg, immediate. = TRUE)
        return(NULL)
      }
      return(model)
    }, error = function(e) {
      msg <- paste("Errore nel modello per il gene:", gene, "-", e$message)
      warning(msg, immediate. = TRUE)
      return(NULL)
    })
  }, warning = function(w) {
    message(w$message)
    invokeRestart("muffleWarning")
  })
}

calculate_mutation_frequency <- function(df, gene_columns) {
  # 1. Calcola le categorie BMI e patient_id
  df <- df %>%
    mutate(
      patient_id = substr(Tumor_Sample_Barcode, 1, 12),
      bmi_sottotipo_4cat = case_when(
        bmi < 18.5 ~ "Underweight",
        bmi >= 18.5 & bmi < 25 ~ "Normal weight",
        bmi >= 25 & bmi < 30 ~ "Overweight",
        bmi >= 30 ~ "Obese",
        TRUE ~ NA_character_
      ),
      bmi_sottotipo = case_when(
        bmi_sottotipo_4cat %in% c("Underweight", "Normal weight") ~ "Normal_Underweight",
        bmi_sottotipo_4cat %in% c("Overweight", "Obese") ~ "Overweight_Obese",
        TRUE ~ NA_character_
      ),
      bmi_sottotipo = factor(bmi_sottotipo,
                             levels = c("Normal_Underweight", "Overweight_Obese"))
    ) %>%
    filter(!is.na(bmi_sottotipo))
  
  # 2. Calcola il numero totale di pazienti per categoria BMI
  total_patients_by_bmi <- df %>%
    distinct(patient_id, bmi_sottotipo) %>%
    count(bmi_sottotipo, name = "total_patients")
  
  # 3. Processa le mutazioni
  result <- df %>%
    select(all_of(gene_columns), patient_id, bmi_sottotipo) %>%
    tidyr::pivot_longer(
      cols = all_of(gene_columns),
      names_to = "gene",
      values_to = "mutation"
    ) %>%
    filter(mutation == 1) %>%
    distinct(gene, patient_id, bmi_sottotipo, .keep_all = TRUE) %>%
    group_by(gene, bmi_sottotipo) %>%
    summarise(
      n_pazienti_con_mutazione = n_distinct(patient_id),
      .groups = "drop"
    ) %>%
    left_join(total_patients_by_bmi, by = "bmi_sottotipo") %>%
    mutate(
      frequenza = n_pazienti_con_mutazione / total_patients,
      sd_frequenza = sqrt(frequenza * (1 - frequenza) / total_patients),
      total_patients = total_patients
    )
  
  return(result)
}


create_mutation_matrix <- function(maf_obj) {
  if (nrow(maf_obj@data) == 0) {
    return(data.frame(Sample = character(), Gene = character(), Mutation = integer()))
  }
  
  mut_matrix <- maf_obj@data %>%
    dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification) %>%
    filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 
                                         'Frame_Shift_Del', 'Frame_Shift_Ins', 
                                         'In_Frame_Del', 'In_Frame_Ins', 
                                         'Splice_Site', 'Translation_Start_Site')) %>%
    distinct(Tumor_Sample_Barcode, Hugo_Symbol, .keep_all = TRUE) %>%
    mutate(Mutation = 1)
  
  mut_matrix <- mut_matrix %>%
    tidyr::pivot_wider(names_from = Hugo_Symbol, values_from = Mutation, 
                       values_fill = list(Mutation = 0)) %>%
    as.data.frame()
  
  return(mut_matrix)
}

prepare_regression_data <- function(df) {
  gene_cols <- names(df)[(which(names(df) == "Variant_Classification") + 1):ncol(df)]
  
  # Crea formula base con variabili disponibili
  base_vars <- "bmi"  # Sempre presente
  if ("age_at_diagnosis" %in% colnames(df) && !all(is.na(df$age_at_diagnosis))) {
    base_vars <- c(base_vars, "age_at_diagnosis")
  }
  if ("gender" %in% colnames(df) && !all(is.na(df$gender))) {
    base_vars <- c(base_vars, "gender")
  }
  if ("tmb" %in% colnames(df) && !all(is.na(df$tmb))) {
    base_vars <- c(base_vars, "tmb")
  }
  
  # Prepara i dati
  df <- df %>%
    mutate(across(all_of(gene_cols), ~factor(., levels = c("0", "1")))) %>%
    mutate(
      Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode),
      Variant_Classification = as.character(Variant_Classification),
      bmi = as.numeric(bmi),
      age_at_diagnosis = as.numeric(age_at_diagnosis),
      gender = as.factor(gender),
      tmb = as.numeric(tmb),
      bcr_patient_barcode = as.character(bcr_patient_barcode)
    )
  
  names(df) <- gsub("-", "_", names(df))
  gene_cols <- gsub("-", "_", gene_cols)
  
  return(list(data = df, gene_cols = gene_cols, base_vars = base_vars))
}

generate_gene_ranking_csv <- function(results_file, frequency_file, output_filename) {
  tryCatch({
    # Leggi i file di input
    regression_results <- read.csv(results_file)
    mutation_frequency <- read.csv(frequency_file)
    
    # Esegui il left_join per la colonna 'gene'
    gene_ranking_df <- left_join(regression_results, mutation_frequency, by = "gene")
    
    # Rinomina alcune colonne per coerenza con la versione precedente (opzionale)
    gene_ranking_df <- gene_ranking_df %>%
      rename(
        num_patients_mutated = n_pazienti_con_mutazione,
        mutation_frequency = frequenza,
        mutation_frequency_sd = sd_frequenza
      )
    
    # Scrivi il dataframe risultante nel file CSV di output
    write.csv(gene_ranking_df, output_filename, row.names = FALSE)
    message("File ", output_filename, " scritto con successo con ", nrow(gene_ranking_df), " righe.")
    
    return(invisible(gene_ranking_df))
    
  }, error = function(e) {
    message("ERRORE in generate_gene_ranking_csv_v2: ", e$message)
    traceback()
    return(NULL)
  })
}

calculate_bmi_categories <- function(clinical_df, output_filename) {
  clinical_df <- clinical_df %>%
    mutate(
      bmi_category_4cat = case_when(
        bmi < 18.5 ~ "Underweight",
        bmi >= 18.5 & bmi < 25 ~ "Normal weight",
        bmi >= 25 & bmi < 30 ~ "Overweight",
        bmi >= 30 ~ "Obese",
        TRUE ~ NA_character_
      ),
      bmi_category = case_when(
        bmi_category_4cat %in% c("Underweight", "Normal weight") ~ "Normal_Underweight",
        bmi_category_4cat %in% c("Overweight", "Obese") ~ "Overweight_Obese",
        TRUE ~ NA_character_
      ),
      bmi_category = factor(bmi_category,
                            levels = c("Normal_Underweight", "Overweight_Obese"))
    ) %>%
    filter(!is.na(bmi_category))
  
  bmi_counts <- clinical_df %>%
    distinct(bcr_patient_barcode, bmi_category) %>%
    count(bmi_category, name = "n_pazienti") %>%
    tidyr::complete(bmi_category, fill = list(n_pazienti = 0))
  
  write.csv(bmi_counts, output_filename, row.names = FALSE)
  
  return(bmi_counts)
}
# Funzione principale -----------------------------------------------------
download_tcga_data <- function(ctype) {
  main_gdc_dir <- "./GDCdata"
  cancer_specific_dir <- file.path(main_gdc_dir, ctype)
  
  # Verifica se la cartella principale GDCdata esiste, altrimenti la crea
  if (!dir.exists(main_gdc_dir)) {
    message("Cartella ", main_gdc_dir, " non trovata. Creazione...")
    dir.create(main_gdc_dir, showWarnings = FALSE)
  }
  
  # Verifica se la cartella specifica del cancro esiste, altrimenti la considera come destinazione del download
  if (!dir.exists(cancer_specific_dir)) {
    message("Cartella specifica per ", ctype, " non trovata in ", main_gdc_dir, ". Avvio del download...")
    tryCatch({
      query <- GDCquery(
        project = ctype,
        data.category = "Simple Nucleotide Variation",
        data.type = "Masked Somatic Mutation",
        workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
      )
      # Scarica direttamente nella cartella specifica del cancro
      GDCdownload(query, method = "api", directory = main_gdc_dir)
      
      # Verifica e sposta la cartella scaricata (se necessario)
      downloaded_folder <- list.files(main_gdc_dir, pattern = paste0("GDCdata_", str_replace_all(ctype, "-", "_")))
      if (length(downloaded_folder) == 1) {
        old_path <- file.path(main_gdc_dir, downloaded_folder)
        new_path <- cancer_specific_dir
        if (file.exists(old_path) && !dir.exists(new_path)) {
          file.rename(old_path, new_path)
          message(paste("Dati TCGA per", ctype, "scaricati e spostati in", new_path))
        } else if (dir.exists(new_path)) {
          message(paste("Cartella di destinazione", new_path, "esiste già."))
          # Potresti voler aggiungere qui una logica per gestire i file esistenti
        } else {
          warning(paste("Impossibile spostare la cartella scaricata per", ctype, ". Percorso sorgente:", old_path, ", percorso destinazione:", new_path))
        }
      } else if (length(downloaded_folder) > 1) {
        warning(paste("Trovate più cartelle di download per", ctype, ". Richiesta revisione."))
      } else {
        warning(paste("Nessuna cartella di download trovata per", ctype, "."))
      }
      
    }, error = function(e) {
      stop(paste("Errore durante il download dei dati TCGA per", ctype, ":", e$message))
    })
  } else {
    message("Cartella specifica per ", ctype, " trovata in ", main_gdc_dir, ". Si procederà con i file presenti in:", cancer_specific_dir)
  }
}
run_analysis <- function(ctype, string) {
  message("Avvio analisi per: ", ctype)
  
  download_tcga_data(ctype)
  # Query e preparazione dati
  TCGAbiolinks:::getProjectSummary(ctype)
  
  query <- GDCquery(
    project = ctype,
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  
  # Prepara dati SNV
  dir_path <- "./GDCdata"
  SNV_data <- GDCprepare(query, directory = "GDCdata", summarizedExperiment = TRUE)
  
  maf_data <- SNV_data
  message(paste("Tentativo di leggere il file maf"))
  maf <- read.maf(maf = maf_data)
  
  # Estrai mutazioni
  nonsynonymous_mut <- subsetMaf(maf, query = "Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 
                                  'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 
                                  'In_Frame_Ins', 'Splice_Site', 'Translation_Start_Site')")
  
  oncogenic_mut <- subsetMaf(maf, query = "Hugo_Symbol %in% driver_genes")
  
  # Carica dati clinici
  directory <- "/media/rb/Data1/lavoro_unipa_stassi/xTania/analisi_TCGA_bmi_x_cnv(Tania)/IN_GDCdata/clinical_data"
  file_path <- list.files(path = directory, pattern = paste0(string, "$"), recursive = TRUE, full.names = TRUE)
  message(paste("Tentativo di leggere il file clinico:", file_path))
  df <- read.table(file_path, header = TRUE, sep = "\t", fill = TRUE)
  parsing_problems <- problems(df)
  if(nrow(parsing_problems) > 0) {
    warning("Problemi di parsing nel file: ", file_path)
    print(parsing_problems)
    # Opzionale: salva i problemi in un file
    write.csv(parsing_problems, paste0(file_path, "_parsing_problems.csv"), row.names = FALSE)
  }
  # Controlla i problemi di parsing subito dopo la lettura
  parsing_problems_clinical <- problems(df)
  if (nrow(parsing_problems_clinical) > 0) {
    warning(paste("Problemi di parsing rilevati nel file clinico:", file_path))
    print(parsing_problems_clinical)
    # Potresti anche decidere di interrompere l'esecuzione qui se i problemi sono gravi:
    # stop("Problemi di parsing critici nel file clinico.")
  }
  df1 <- df[-c(1, 2), ]
  
  # Gestione colonne mancanti
  required_columns <- c("bcr_patient_barcode", "weight_kg_at_diagnosis", "height_cm_at_diagnosis")
  optional_columns <- c("age_at_diagnosis", "gender")
  
  missing_required <- setdiff(required_columns, colnames(df1))
  if (length(missing_required) > 0) {
    stop(paste("Colonne obbligatorie mancanti:", paste(missing_required, collapse=", ")))
  }
  
  for (col in optional_columns) {
    if (!col %in% colnames(df1)) {
      df1[[col]] <- NA
      message("Colonna opzionale ", col, " non trovata - impostata a NA per ", string)
    }
  }
  
  # Calcola BMI
  df1 <- df1 %>%
    mutate(
      weight_kg_at_diagnosis = as.numeric(gsub("[^0-9.]", NA, weight_kg_at_diagnosis)),
      height_cm_at_diagnosis = as.numeric(gsub("[^0-9.]", NA, height_cm_at_diagnosis)),
      bmi = ifelse(height_cm_at_diagnosis > 0 & weight_kg_at_diagnosis > 0,
                   weight_kg_at_diagnosis / ((height_cm_at_diagnosis/100)^2),
                   NA_real_),
      age_at_diagnosis = as.numeric(ifelse(grepl("^[0-9.]+$", age_at_diagnosis),
                                           age_at_diagnosis,
                                           NA_real_))
    )
  
  df_senza_na <- df1 %>% filter(!is.na(bmi))
  clinical_df <- df_senza_na
  # **INSERISCI QUI IL CALCOLO E IL SALVATAGGIO DELLE CATEGORIE BMI**
  calculate_bmi_categories(clinical_df, paste0(string, "_bmi_categories.csv"))
  message("Calcolo categorie BMI completato per: ", ctype)
  
  # Calcola TMB
  callable_genome_size_mb <- 40
  tmb_results <- tmb(maf, captureSize = callable_genome_size_mb)
  tmb_results$Sample_cut <- substr(tmb_results$Tumor_Sample_Barcode, 1, 12)
  
  clinical_df <- merge(clinical_df, tmb_results[, c("Sample_cut", "total")], 
                       by.x = "bcr_patient_barcode", 
                       by.y = "Sample_cut", 
                       all.x = TRUE) %>%
    dplyr::rename(tmb = total) %>%
    mutate(tmb = ifelse(is.na(tmb), 0, tmb))
  
  # Crea matrici mutazioni
  mutation_matrix_nonsynonymous <- create_mutation_matrix(nonsynonymous_mut)
  mutation_matrix_oncogenic <- create_mutation_matrix(oncogenic_mut)
  
  # Merge con dati clinici
  mutation_matrix_nonsynonymous$bcr_patient_barcode <- substr(mutation_matrix_nonsynonymous$Tumor_Sample_Barcode, 1, 12)
  merged_df_nonsynonymous <- merge(clinical_df, mutation_matrix_nonsynonymous, 
                                   by.x = "bcr_patient_barcode", 
                                   by.y = "bcr_patient_barcode") %>%
    replace(is.na(.), 0)
  
  mutation_matrix_oncogenic$bcr_patient_barcode <- substr(mutation_matrix_oncogenic$Tumor_Sample_Barcode, 1, 12)
  merged_df_oncogenic <- merge(clinical_df, mutation_matrix_oncogenic, 
                               by.x = "bcr_patient_barcode", 
                               by.y = "bcr_patient_barcode") %>%
    replace(is.na(.), 0)
  
  # Filtra per BMI
  bmi_min <- 10
  bmi_max <- 60
  merged_df_nonsynonymous <- merged_df_nonsynonymous %>% filter(bmi >= bmi_min, bmi <= bmi_max)
  merged_df_oncogenic <- merged_df_oncogenic %>% filter(bmi >= bmi_min, bmi <= bmi_max)
  
  # Calcola frequenze mutazioni
  gene_columns_ns <- names(merged_df_nonsynonymous)[(which(names(merged_df_nonsynonymous) == "Variant_Classification") + 1):ncol(merged_df_nonsynonymous)]
  mutation_frequency_ns <- calculate_mutation_frequency(merged_df_nonsynonymous, gene_columns_ns)
  write.csv(mutation_frequency_ns, paste0(string, "_mutation_frequency_ns_freq.csv"), row.names = FALSE)
  
  gene_columns_onc <- names(merged_df_oncogenic)[(which(names(merged_df_oncogenic) == "Variant_Classification") + 1):ncol(merged_df_oncogenic)]
  mutation_frequency_onc <- calculate_mutation_frequency(merged_df_oncogenic, gene_columns_onc)
  write.csv(mutation_frequency_onc, paste0(string, "_mutation_frequency_onc_freq.csv"), row.names = FALSE)
  
  # Prepara dati per regressione nonsynonymous
  prep_ns <- prepare_regression_data(merged_df_nonsynonymous)
  merged_df_ns <- prep_ns$data
  gene_cols_ns <- prep_ns$gene_cols
  
  # Regressione semplice (solo BMI) per nonsynonymous
  results_ns <- data.frame()
  for (gene in gene_cols_ns) {
    model_ns <- run_regression(merged_df_ns, gene, "bmi")
    
    if (!is.null(model_ns)) {
      summary_ns <- summary(model_ns)
      if ("bmi" %in% rownames(summary_ns$coefficients)) {
        results_ns <- rbind(results_ns, data.frame(
          gene = gene,
          estimate = summary_ns$coefficients["bmi", "Estimate"],
          p.value = summary_ns$coefficients["bmi", "Pr(>|z|)"],
          status = ifelse(model_ns$converged, "success", "failed")
        ))
      }
    }
  }
  
  results_ns$p.adj <- p.adjust(results_ns$p.value, method = "BH")
  results_ns_filtered <- results_ns %>% filter(p.adj != 0)
  write.csv(results_ns_filtered, paste0(string, "_results_ns.csv"), row.names = FALSE)
  if (nrow(results_ns_filtered) > 0) {
    generate_gene_ranking_csv(
      results_file = paste0(string, "_results_ns.csv"),
      frequency_file = paste0(string, "_mutation_frequency_ns_freq.csv"),
      output_filename = paste0(string, "_gene_ranking_ns.csv")
    )
  }
  
  
  # Regressione completa per nonsynonymous
  results_ns_full <- data.frame()
  if (length(prep_ns$base_vars) > 1) {
    for (gene in prep_ns$gene_cols) {
      model_ns_full <- run_regression(merged_df_ns, gene, prep_ns$base_vars)
      
      if (!is.null(model_ns_full)) {
        coef_summary <- summary(model_ns_full)$coefficients
        temp_results <- data.frame(gene = gene)
        
        for (var in prep_ns$base_vars) {
          if (var %in% rownames(coef_summary)) {
            temp_results[[paste0("estimate_", var)]] <- coef_summary[var, "Estimate"]
            temp_results[[paste0("p.value_", var)]] <- coef_summary[var, "Pr(>|z|)"]
          }
        }
        results_ns_full <- bind_rows(results_ns_full, temp_results)
      }
    }
    
    for (var in prep_ns$base_vars) {
      p_col <- paste0("p.value_", var)
      padj_col <- paste0("padj_", var)
      
      if (p_col %in% names(results_ns_full)) {
        results_ns_full[[padj_col]] <- p.adjust(results_ns_full[[p_col]], method = "BH")
      }
    }
  } else {
    message(paste(string, "- Nessuna variabile opzionale disponibile per la regressione completa (nonsynonymous)."))
  }
  
  if (nrow(results_ns_full) > 0) {
    write.csv(results_ns_full, paste0(string, "_results_ns_full.csv"), row.names = FALSE)
  } else {
    message("Nessun risultato per la regressione completa (nonsynonymous) in ", string)
    write.csv(data.frame(), paste0(string, "_results_ns_full.csv"), row.names = FALSE)
  }
  
  # Ripeti lo stesso processo per le mutazioni oncogeniche
  prep_onc <- prepare_regression_data(merged_df_oncogenic)
  merged_df_onc <- prep_onc$data
  gene_cols_onc <- prep_onc$gene_cols
  
  # Regressione semplice (solo BMI) per oncogeniche
  results_onc <- data.frame()
  for (gene in gene_cols_onc) {
    model_onc <- run_regression(merged_df_onc, gene, "bmi")
    
    if (!is.null(model_onc)) {
      summary_onc <- summary(model_onc)
      if ("bmi" %in% rownames(summary_onc$coefficients)) {
        results_onc <- rbind(results_onc, data.frame(
          gene = gene,
          estimate = summary_onc$coefficients["bmi", "Estimate"],
          p.value = summary_onc$coefficients["bmi", "Pr(>|z|)"],
          status = ifelse(model_onc$converged, "success", "failed")
        ))
      }
    }
  }
  
  results_onc$p.adj <- p.adjust(results_onc$p.value, method = "BH")
  results_onc_filtered <- results_onc %>% filter(p.adj != 0)
  write.csv(results_onc_filtered, paste0(string, "_results_onc.csv"), row.names = FALSE)
  if (nrow(results_onc_filtered) > 0) {
    generate_gene_ranking_csv(
      results_file = paste0(string, "_results_onc.csv"),
      frequency_file = paste0(string, "_mutation_frequency_onc_freq.csv"),
      output_filename = paste0(string, "_gene_ranking_onc.csv")
    )
  }
  
  
  # Regressione completa per oncogeniche
  # Regressione completa per oncogeniche
  results_onc_full <- data.frame()
  if (length(prep_onc$base_vars) > 1) {
    for (gene in prep_onc$gene_cols) {
      model_onc_full <- run_regression(merged_df_onc, gene, prep_onc$base_vars)
      
      if (!is.null(model_onc_full)) {
        coef_summary <- summary(model_onc_full)$coefficients
        temp_results <- data.frame(gene = gene)
        
        for (var in prep_onc$base_vars) {
          if (var %in% rownames(coef_summary)) {
            temp_results[[paste0("estimate_", var)]] <- coef_summary[var, "Estimate"]
            temp_results[[paste0("p.value_", var)]] <- coef_summary[var, "Pr(>|z|)"]
          }
        }
        results_onc_full <- bind_rows(results_onc_full, temp_results)
      }
    }
    
    for (var in prep_onc$base_vars) {
      p_col <- paste0("p.value_", var)
      padj_col <- paste0("padj_", var)
      
      if (p_col %in% names(results_onc_full)) {
        results_onc_full[[padj_col]] <- p.adjust(results_onc_full[[p_col]], method = "BH")
      }
    }
  } else {
    message(paste(string, "- Nessuna variabile opzionale disponibile per la regressione completa (oncogeniche)."))
  }
  
  if (nrow(results_onc_full) > 0) {
    write.csv(results_onc_full, paste0(string, "_results_onc_full.csv"), row.names = FALSE)
  } else {
    message("Nessun risultato per la regressione completa (oncogeniche) in ", string)
    write.csv(data.frame(), paste0(string, "_results_onc_full.csv"), row.names = FALSE)
  }
  
  # Visualizzazione TMB
  winsorize_tmb <- function(x, probs = c(0.05, 0.95)) {
    q <- quantile(x, probs, na.rm = TRUE)
    x[x < q[1]] <- q[1]
    x[x > q[2]] <- q[2]
    return(x)
  }
  
  clinical_df$tmb_winsorized <- winsorize_tmb(clinical_df$tmb)
  
  # Define the new BMI categories with desired labels
  clinical_df$bmi_category_4cat <- case_when(
    clinical_df$bmi < 18.5 ~ "Underweight",
    clinical_df$bmi >= 18.5 & clinical_df$bmi < 25 ~ "Normal weight",
    clinical_df$bmi >= 25 & clinical_df$bmi < 30 ~ "Overweight",
    clinical_df$bmi >= 30 ~ "Obese",
    TRUE ~ NA_character_
  )
  clinical_df$bmi_category <- factor(
    case_when(
      clinical_df$bmi_category_4cat %in% c("Underweight", "Normal weight") ~ "Normal_Underweight",
      clinical_df$bmi_category_4cat %in% c("Overweight", "Obese") ~ "Overweight_Obese",
      TRUE ~ NA_character_
    ),
    levels = c("Normal_Underweight", "Overweight_Obese"),
    labels = c("BMI < 25", "BMI \u2265 25") # New labels here!
  )
  
  # Plot TMB (Winsorized) with standardized theme and colors
  plot_tmb_bmi <- ggplot(clinical_df, aes(x = bmi_category, y = tmb_winsorized, fill = bmi_category, color = bmi_category)) +
    geom_violin(alpha = 0.6, linewidth = 1) + # Increased linewidth for better visibility
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7) + # Slightly increased jitter point size
    geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) + # Added outline to boxplot, removed outliers
    labs(title = paste("Winsorized TMB Distribution for", string), # Adjusted title
         x = "BMI Category",
         y = "Winsorized TMB") +
    scale_fill_manual(values = c("BMI < 25" = "orange3", "BMI \u2265 25" = "dodgerblue4")) +
    scale_color_manual(values = c("BMI < 25" = "orange3", "BMI \u2265 25" = "dodgerblue4")) +
    theme_minimal() + # Start with minimal theme for full control
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 64, family = "Arial", color = "black"),
      axis.title.x = element_text(size = 64, family = "Arial", color = "black", face = "bold", margin = margin(t = 20)),
      axis.title.y = element_text(size = 64, family = "Arial", color = "black", face = "bold", margin = margin(r = 20)),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 57, family = "Arial", color = "black", face = "bold"),
      axis.text.y = element_text(size = 57, family = "Arial", color = "black", face = "bold"),
      legend.position = "none", # Removed legend as colors are explicit for categories
      panel.grid.major = element_blank(), # No major grid
      panel.grid.minor = element_blank(), # No minor grid
      panel.background = element_rect(fill = "white", color = NA), # White background
      plot.background = element_rect(fill = "white", color = NA), # White plot background
      axis.line = element_line(color = "black") # Black axis lines
    )
  
  ggsave(paste0(string, "_tmb_bmi_winsorized_violin_standardized.tiff"), plot_tmb_bmi, width = 10, height = 8, dpi = 300, bg = "white") # Save as TIFF
  
  # Plot TMB (No Winsorized) with standardized theme and colors
  plot_tmb_bmi_no_winsor <- ggplot(clinical_df, aes(x = bmi_category, y = tmb, fill = bmi_category, color = bmi_category)) +
    geom_violin(alpha = 0.6, linewidth = 1) + # Increased linewidth for better visibility
    geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.7) + # Slightly increased jitter point size
    geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) + # Added outline to boxplot, removed outliers
    labs(title = paste("TMB Distribution for", string), # Adjusted title
         x = "BMI Category",
         y = "TMB") +
    scale_fill_manual(values = c("BMI < 25" = "orange3", "BMI \u2265 25" = "dodgerblue4")) +
    scale_color_manual(values = c("BMI < 25" = "orange3", "BMI \u2265 25" = "dodgerblue4")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 64, family = "Arial", color = "black"),
      axis.title.x = element_text(size = 64, family = "Arial", color = "black", face = "bold", margin = margin(t = 20)),
      axis.title.y = element_text(size = 64, family = "Arial", color = "black", face = "bold", margin = margin(r = 20)),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 57, family = "Arial", color = "black", face = "bold"),
      axis.text.y = element_text(size = 57, family = "Arial", color = "black", face = "bold"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black")
    )
  
  ggsave(paste0(string, "_tmb_bmi_violin_no_winsor_standardized.tiff"), plot_tmb_bmi_no_winsor, width = 10, height = 8, dpi = 300, bg = "white") # Save as TIFF
}

# Esecuzione --------------------------------------------------------------
input_file <- "cancer_types.tsv"
message(paste("Tentativo di leggere il file cancer_types:", input_file))
cancer_data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Note: 'problems(df)' here would refer to a 'df' variable which is not defined in this scope.
# The 'problems' function is part of readr and typically used right after a read_csv/tsv call.
# It's okay to remove or comment this specific check if it's not applicable here.
# parsing_problems <- problems(df) 
# if(nrow(parsing_problems) > 0) {
#   warning("Problemi di parsing nel file: ", file_path)
#   print(parsing_problems)
#   write.csv(parsing_problems, paste0(file_path, "_parsing_problems.csv"), row.names = FALSE)
# }

for (i in 1:nrow(cancer_data)) {
  ctype <- cancer_data[i, "ctype"]
  string <- cancer_data[i, "string"]
  
  message("\n=== Inizio analisi per: ", ctype, " ===")
  
  tryCatch({
    run_analysis(ctype, string)
    message("Analisi completata con successo per: ", ctype)
  }, error = function(e) {
    message("ERRORE in ", ctype, ": ", e$message)
  })
}

# Chiudi i sink
# Chiudi i sink in modo robusto alla fine dell'analisi
tryCatch({
  # Chiudi sink standard
  attempts <- 0
  while (sink.number() > 0 && attempts < 10) {
    sink()
    attempts <- attempts + 1
    Sys.sleep(0.1)
  }
  
  # Inserisci qui il controllo per chiudere warn_con
  if (exists("warn_con") && isOpen(warn_con)) {
    try(close(warn_con), silent = TRUE)
  }
  
  # Chiudi sink message
  attempts <- 0
  while (sink.number(type = "message") > 0 && attempts < 10) {
    try(sink(type = "message"), silent = TRUE)
    attempts <- attempts + 1
    Sys.sleep(0.1)
  }
  
  # Verifica
  if (sink.number() > 0 || sink.number(type = "message") > 0) {
    warning("Alcuni sink potrebbero non essere stati chiusi correttamente")
  }
  
  # Stampa eventuali warning
  if (file.exists(warn_file)) {
    warn_lines <- readLines(warn_file)
    if (length(warn_lines) > 0) {
      message("--- WARNINGS ---")
      message(paste(warn_lines, collapse = "\n"))
    }
  }
}, error = function(e) {
  message("Errore durante la chiusura finale dei sink: ", e$message)
})