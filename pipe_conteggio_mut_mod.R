library(TCGAbiolinks)
library(maftools)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(viridis)
library(scales)
library(patchwork)
library(RColorBrewer)
library(purrr)
library(randomcoloR)
library(showtext) # Aggiunto per gestione font

# Attiva il rendering dei font di sistema con showtext
showtext_auto()

# Funzione per scurire i colori (mantenuta come nel tuo script originale)
darken_rgb <- function(color, factor = 0.5) {
  col <- col2rgb(color)
  alpha <- alpha(color)
  if (is.na(alpha)) {
    new_col <- pmax(col * factor, 0)
    return(rgb(new_col[1, ], new_col[2, ], new_col[3, ], maxColorValue = 255))
  } else {
    new_col <- pmax(col * factor, 0)
    return(rgb(new_col[1, ], new_col[2, ], new_col[3, ], alpha = alpha * 255, maxColorValue = 255))
  }
}

# 1. Definisci il tipo di cancro (BLCA)
ctype <- "TCGA-BLCA"

# 2. Definisci il percorso della cartella GDCdata
dir_path <- "./GDCdata"

# 3. Carica i dati di mutazione
query_mut <- GDCquery(
  project = ctype,
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
SNV_data <- GDCprepare(query_mut, directory = dir_path, summarizedExperiment = FALSE)
BLCA_maf <- read.maf(SNV_data)

# 4. Carica il file clinico
directory_clinico <- "/media/rb/Data1/lavoro_unipa_stassi/xTania/analisi_TCGA_bmi_x_cnv(Tania)/IN_GDCdata/clinical_data/gdc_download_20250131_134224.911560/eaa71705-960a-4abd-b5d7-f5fdc0d0c5af/"
string_file_clinico <- "nationwidechildrens.org_clinical_patient_blca.txt"

# CORREZIONE: Utilizza list.files per ottenere il percorso del file come stringa
file_path_clinico <- list.files(path = directory_clinico, pattern = paste0(string_file_clinico, "$"), recursive = TRUE, full.names = TRUE)
clinical_df <- read.table(file_path_clinico, header = TRUE, sep = "\t", fill = TRUE)
parsing_problems_clinical <- problems(clinical_df)
if (nrow(parsing_problems_clinical) > 0) print(parsing_problems_clinical)
clinical_df <- clinical_df[-c(1, 2), ]
clinical_df <- clinical_df %>% rename(Tumor_Sample_Barcode = bcr_patient_barcode)

# 5. Calcola il BMI
clinical_bmi <- clinical_df %>%
  select(Tumor_Sample_Barcode, weight_kg_at_diagnosis, height_cm_at_diagnosis) %>%
  filter(!is.na(weight_kg_at_diagnosis) & !is.na(height_cm_at_diagnosis)) %>%
  mutate(
    weight_kg_at_diagnosis = as.numeric(gsub("[^0-9.]", NA, weight_kg_at_diagnosis)),
    height_cm_at_diagnosis = as.numeric(gsub("[^0-9.]", NA, height_cm_at_diagnosis)),
    BMI = weight_kg_at_diagnosis / ((height_cm_at_diagnosis/100)^2)) %>%
  filter(!is.na(BMI)) %>%
  mutate(bmi_category = case_when(BMI < 25 ~ "Normal_Under", BMI >= 25 ~ "Over_Obese", TRUE ~ NA_character_)) %>%
  select(Tumor_Sample_Barcode, BMI, bmi_category)

# 6. Trunca Tumor_Sample_Barcode in BLCA_maf@data
BLCA_maf@data <- BLCA_maf@data %>%
  mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 12))

# Crea la variabile lista_completa_geni
# ORDINE DEI GENI AGGIORNATO
lista_completa_geni <- c("BRCA2", "DNAH9", "GRIA4", "PLXNA4","UNC13C", "FCGBP", "SF3B1","ELP1", "NES", "TRERF1")

# Crea un dataframe vuoto per salvare i risultati di tutti i geni
all_gene_mutation_frequency <- data.frame()

# Itera attraverso la lista dei geni
for (gene in lista_completa_geni) {
  # 7. Estrai le mutazioni non sinonime del gene corrente e unisci con le categorie BMI
  gene_mutations <- BLCA_maf@data %>%
    filter(Hugo_Symbol == gene, Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site")) %>%
    left_join(clinical_bmi, by = "Tumor_Sample_Barcode")
  
  # 8. Estrai l'informazione sull'esone dal MAF object per il gene corrente
  exon_info_gene <- BLCA_maf@data %>%
    filter(Hugo_Symbol == gene, Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site")) %>%
    select(HGVSc, EXON) %>%
    distinct()
  
  # 9. Unisci l'informazione sull'esone al dataframe delle mutazioni del gene corrente
  gene_mutations_with_exon <- gene_mutations %>%
    left_join(exon_info_gene, by = "HGVSc")
  
  # DEDUPLICAZIONE AGGRESSIVA BASATA SULLE COLONNE DI RAGGRUPPAMENTO
  gene_mutations_with_exon_deduped <- gene_mutations_with_exon %>%
    distinct(Tumor_Sample_Barcode, HGVSc, EXON.x, .keep_all = TRUE)
  
  # 1. Calcola la frequenza delle mutazioni per il gene corrente
  mutation_frequency_gene <- gene_mutations_with_exon_deduped %>%
    group_by(bmi_category, EXON.x) %>% # Raggruppa per esone e BMI
    summarise(
      patients_with_mutations = n(), # Conta i pazienti con mutazioni
      patients = paste(unique(Tumor_Sample_Barcode), collapse = ", ")
    ) %>%
    arrange(desc(patients_with_mutations)) %>%
    mutate(gene = gene, EXON = EXON.x)
  
  # Aggiungi i risultati del gene corrente al dataframe di tutti i geni
  all_gene_mutation_frequency <- bind_rows(all_gene_mutation_frequency, mutation_frequency_gene)
}

# Rinomina la colonna EXON.x e rimuovi EXON.x
all_gene_mutation_frequency <- all_gene_mutation_frequency %>%
  select(-EXON.x)

# 11. Calcola il numero totale di pazienti per categoria BMI
total_patients_by_bmi <- clinical_bmi %>%
  group_by(bmi_category) %>%
  summarise(total_patients = n())

# 12. Unisci il numero totale di pazienti per BMI al dataframe delle frequenze
all_gene_mutation_frequency_with_totals <- all_gene_mutation_frequency %>%
  left_join(total_patients_by_bmi, by = "bmi_category")

# 13. Normalizza il conteggio delle mutazioni
all_gene_mutation_frequency_normalized <- all_gene_mutation_frequency_with_totals %>%
  mutate(
    normalized_mutation_frequency = patients_with_mutations / total_patients
  )

# 14. Raggruppa gli esoni per visualizzazione
all_gene_mutation_frequency_grouped <- all_gene_mutation_frequency_normalized %>%
  group_by(gene) %>%
  mutate(total_exon_count = sum(patients_with_mutations, na.rm = TRUE)) %>% # Calcola il conteggio totale per esone
  ungroup() %>%
  mutate(EXON_DISPLAY = ifelse(total_exon_count == 1, "Other Exons", as.character(EXON))) %>%
  group_by(gene, bmi_category, EXON_DISPLAY) %>%
  summarise(
    normalized_mutation_frequency = mean(normalized_mutation_frequency, na.rm = TRUE), # Usa la media della frequenza normalizzata
    patients = paste(unique(unlist(strsplit(patients, ", "))), collapse = ", ")
  ) %>%
  ungroup()

# FILTRAGGIO PER RIMUOVERE RIGHE CON NA IN QUALSIER COLONNA
all_gene_mutation_frequency_grouped <- all_gene_mutation_frequency_grouped %>%
  filter(complete.cases(.))

# Preparazione dei dati per il plot
df_plot <- all_gene_mutation_frequency_grouped %>%
  filter(
    !is.na(EXON_DISPLAY),
    bmi_category %in% c("Over_Obese", "Normal_Under")
  ) %>%
  mutate(
    bmi_category = case_when(
      bmi_category == "Normal_Under" ~ "<25",
      bmi_category == "Over_Obese" ~ "≥25",
      TRUE ~ bmi_category
    ),
    gene_exon_label = paste0(gene, "_", EXON_DISPLAY)
  )

# Rimuovi il testo a partire dal simbolo '/' dalle etichette degli esoni
df_plot$EXON_DISPLAY <- sub("/.*", "", df_plot$EXON_DISPLAY)


# Calcola l'ordine degli esoni a livello globale, gestendo i NA
exon_order <- df_plot %>%
  group_by(gene_exon_label) %>%
  summarise(total_frequency = mean(normalized_mutation_frequency, na.rm = TRUE)) %>% # Usa la media della frequenza normalizzata
  arrange(desc(total_frequency)) %>%
  pull(gene_exon_label) %>%
  unique()

# Definisci l'ordine dei geni esattamente come richiesto
gene_order_for_plot <- c("BRCA2", "DNAH9", "GRIA4", "PLXNA4","UNC13C", "FCGBP", "SF3B1","ELP1", "NES", "TRERF1")

# Applica l'ordine dei geni in df_plot (invertito per visualizzazione da alto a basso in orizzontale)
df_plot <- df_plot %>%
  mutate(
    gene = factor(gene, levels = rev(gene_order_for_plot)) # Invertito l'ordine per visualizzazione orizzontale da alto a basso
  )

# Calcola il numero di esoni per ogni gene (necessario per la colorazione)
gene_exon_counts <- df_plot %>%
  group_by(gene) %>%
  summarise(exon_count = n_distinct(EXON_DISPLAY)) %>%
  ungroup()

# Mappiamo i colori per ogni esone (ogni gene ha colori distinti per i suoi esoni)
n_distinct_exons <- length(unique(df_plot$gene_exon_label))
custom_palette <- distinctColorPalette(n_distinct_exons)
gene_exon_color_mapping <- data.frame(
  gene_exon_label = unique(df_plot$gene_exon_label),
  gene_exon_color = custom_palette[1:n_distinct_exons]
)

# Uniamo la mappatura dei colori agli esoni al dataframe df_plot
df_plot <- df_plot %>%
  left_join(gene_exon_color_mapping, by = "gene_exon_label")

# Converti in lista per ggplot
exon_color_list <- setNames(df_plot$gene_exon_color, df_plot$gene_exon_label)

# Applica l'ordine degli esoni per la colorazione (non è influenzato dalla rotazione dell'asse)
df_plot <- df_plot %>%
  mutate(
    gene_exon_label = factor(gene_exon_label, levels = exon_order)
  )

# Crea i grafici
plot_list <- list()
bmi_categories <- unique(df_plot$bmi_category)
y_limit <- 0.15 # Definisci il limite superiore per l'asse y

for (bmi_cat in bmi_categories) {
  df_subset <- df_plot %>% filter(bmi_category == bmi_cat)
  
  p <- ggplot(df_subset, aes(x = gene, y = normalized_mutation_frequency, fill = gene_exon_label)) + # X è il gene, Y è la frequenza
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.5, width = 0.9) + # Larghezza delle barre mantenuta a 0.5
    # Aggiungi le etichette degli esoni direttamente sulle barre
    geom_text(
      aes(label = EXON_DISPLAY), # Solo il nome dell'esone
      position = position_stack(vjust = 0.5), # Centra il testo all'interno di ogni segmento
      color = "black", # Colore del testo
      size = 22, # Dimensione del testo aumentata di 2 punti (da 20 a 22)
      family = "Arial", # Imposta il font Arial per le etichette degli esoni
      fontface = "bold", # Imposta il testo in grassetto per le etichette degli esoni
      # Filtra per mostrare etichette solo su segmenti sufficientemente grandi
      data = df_subset %>% filter(normalized_mutation_frequency > 0.002) # Soglia abbassata per mostrare più etichette
    ) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01), limits = c(0, y_limit)) + # Scala per l'asse y (numerico)
    scale_fill_manual(values = exon_color_list) + # Usa la lista dei colori per la mappatura
    theme_bw() +
    theme(
      # Impostazioni del tema per allineare al grafico precedente
      axis.text.x = element_text(size = 77, family = "Arial", face = "bold"),
      axis.text.y = element_text(size = 77, family = "Arial", face = "bold", angle = 0, hjust = 1),
      axis.title.x = element_text(size = 82, face = "bold", family = "Arial"),
      axis.title.y = element_text(size = 82, face = "bold", family = "Arial"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 82, family = "Arial"), # Titolo del singolo BMI panel
      strip.text = element_text(face = "bold", size = 62, family = "Arial"), # Adjust strip text size if needed (e.g., if you add facets here)
      panel.grid.major = element_blank(), # Rimuovi la griglia maggiore
      panel.grid.minor = element_blank(), # Rimuovi la griglia minore
      legend.position = "none" # Rimuovi la legenda
    ) +
    coord_flip() # Ruota il grafico di 90 gradi per renderlo orizzontale
  
  plot_list[[bmi_cat]] <- p
}

# Combinazione dei grafici usando patchwork
combined_plot <- plot_list[["<25"]] / plot_list[["≥25"]] +
  plot_layout(guides = 'collect')

# Aggiungi il titolo principale separatamente
combined_plot <- combined_plot + plot_annotation(
#  title = "Normalized Mutation Frequency by Gene and Exon in BLCA Patients Stratified by BMI Category",
  theme = theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 72, family = "Arial") # Overall plot title
  )
) &
  theme(
    legend.position = "none",
    # These global theme settings apply to all combined plots.
    # We already set face and family for axis.title/text in the loop,
    # but repeating ensures consistency or acts as fallback.
    axis.title = element_text(face = "bold", family = "Arial"),
    axis.text = element_text(face = "bold", family = "Arial"),
    strip.text = element_text(face = "bold", family = "Arial") # For potential future facets within sub-plots
  )

# Salva il grafico in formato TIFF
ggsave(
  filename = "normalized_mutation_frequency_by_gene_exon_standardized.tiff", # Updated filename
  plot = combined_plot,
  width = 15, # Consistent with previous plots (e.g., heatmap)
  height = 18, # Adjusted for desired elongated shape
  dpi = 300,
  device = "tiff",
  bg = "white"
)

message(paste("Grafico salvato con successo in formato TIFF con scala fissa sull'asse y (0 a", y_limit, "), geni ordinati sull'asse x, etichette degli esoni all'interno delle barre, font Arial, e senza griglia di sfondo. Formattazione allineata al grafico precedente, con titolo principale ora in grassetto."))