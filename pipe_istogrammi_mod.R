library(ggplot2)
library(stringr)
library(dplyr)
library(showtext) # Per gestire i font
library(ggtext)    # Per rendering di testo Markdown/HTML

# Attiva il rendering dei font di sistema con showtext.
showtext_auto()

# Specifica il tipo di tumore
tipo_cancro_specifico <- "BLCA"

# Specifica la lista di geni di interesse
geni_di_interesse <- c("BRCA2", "DNAH9", "GRIA4", "PLXNA4", "UNC13C", "FCGBP", "SF3B1", "ELP1", "NES", "TRERF1")

# Funzione modificata per creare l'unico grafico affiancato (corretto)
crea_unico_grafico_affiancato_corretto <- function(df2_file, geni_lista, tipo_cancro, output_file = "mutation_frequency_genes_bmi_single.tiff") {
  # Carica il dataframe df2
  df2 <- read.csv(df2_file)
  
  # Mappa i nomi BMI in df2 all'inglese
  bmi_mapping <- c(
    "Normal_Underweight" = "Normal_Underweight",
    "Overweight_Obese"   = "Overweight_Obese"
  )
  df2$bmi_sottotipo_en <- factor(
    bmi_mapping[as.character(df2$bmi_sottotipo)],
    levels = bmi_mapping[c("Normal_Underweight", "Overweight_Obese")]
  )
  
  # Filtra df2 per i geni di interesse
  df2_filtrato <- df2 %>% filter(gene %in% geni_lista)
  
  # Sostituisci NA con 0 per frequenza e sd_frequenza
  df2_filtrato$frequenza <- ifelse(is.na(df2_filtrato$frequenza), 0, df2_filtrato$frequenza)
  df2_filtrato$sd_frequenza <- ifelse(is.na(df2_filtrato$sd_frequenza), 0, df2_filtrato$sd_frequenza)
  
  # Ordina i geni come nella lista 'geni_di_interesse'
  df2_filtrato$gene <- factor(df2_filtrato$gene, levels = geni_lista)
  
  # Crea l'unico grafico affiancato (corretto)
  p <- ggplot(df2_filtrato, aes(x = gene, y = frequenza, fill = bmi_sottotipo_en)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = frequenza - sd_frequenza, ymax = frequenza + sd_frequenza),
                  width = 0.2, position = position_dodge(width = 0.9)) +
    scale_fill_manual(
      values = c("Normal_Underweight" = "orange3", "Overweight_Obese" = "dodgerblue4"),
      labels = c("Normal_Underweight" = "<b>< 25</b>", "Overweight_Obese" = "<b>≥ 25</b>"),
      name = "BMI Subtype"
    ) +
    labs(
      x = "Gene",
      y = "Mutation Frequency"
    ) +
    theme_minimal() +
    theme(
      # Standardized font sizes to match volcano plot
      axis.title = element_text(size = 64, family = "Arial", color = "black", face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 69, family = "Arial", color = "black", face = "bold"),
      axis.text.y = element_text(size = 69, family = "Arial", color = "black", face = "bold"),
      legend.title = element_text(size = 57, family = "Arial", color = "black", face = "bold"),
      legend.text = element_markdown(size = 57, family = "Arial", color = "black"),
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black")
    )
  
  # Salva l'unico grafico con dimensioni standardizzate
  ggsave(output_file, plot = p, width = 10, height = 8, units = "in", bg = "white", dpi = 300)
  message(paste("Grafico affiancato unico (corretto) creato per", tipo_cancro, "nonsynonymous e i geni specificati in:", output_file))
}

# Trova e processa i file
df1_files <- list.files(pattern = ".txt_results_ns\\.csv$", full.names = TRUE)

for (df1_file in df1_files) {
  # Estrai il nome del tumore dal nome del file df1
  file_name_base_df1 <- gsub("\\.txt_results_ns\\.csv$", "", basename(df1_file))
  
  if (grepl(tipo_cancro_specifico, file_name_base_df1, ignore.case = TRUE)) {
    # Costruisci il nome del file df2 corrispondente (solo nonsynonymous)
    df2_file_pattern <- paste0(file_name_base_df1, "\\.txt_mutation_frequency_ns_freq\\.csv$")
    df2_files_found <- list.files(path = dirname(df1_file), pattern = df2_file_pattern, full.names = TRUE)
    
    # Debug
    print(paste("df1:", basename(df1_file), "tipo_cancro_file:", file_name_base_df1, "Pattern df2:", df2_file_pattern))
    
    if (length(df2_files_found) == 1) {
      # Chiama la funzione corretta per il grafico unico affiancato per BMI all'interno di ogni gene
      crea_unico_grafico_affiancato_corretto(df2_files_found[1], geni_di_interesse, tipo_cancro_specifico,
                                             output_file = paste0("mutation_frequency_", tipo_cancro_specifico, "_bar_plot.tiff"))
    } else if (length(df2_files_found) > 1) {
      message(paste("Trovati più file df2 nonsynonymous per", tipo_cancro_specifico, ":", paste(basename(df2_files_found), collapse = ", ")))
    } else {
      message(paste("File df2 nonsynonymous non trovato per", tipo_cancro_specifico, "il pattern cercato era:", df2_file_pattern))
    }
  } else {
    message(paste("File df1", basename(df1_file), "non relativo al tumore", tipo_cancro_specifico, ". Ignorato."))
  }
}