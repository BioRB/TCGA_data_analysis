library(ggplot2)
library(dplyr)
library(tidyr)
library(showtext) # Per gestire i font
library(readxl)    # Per leggere file .xls/.xlsx

# Attiva il rendering dei font di sistema con showtext
showtext_auto()

# Definisci la lista di geni di interesse (usata anche qui per consistenza nell'ordine)
geni_di_interesse <- c("BRCA2", "DNAH9", "GRIA4", "PLXNA4", "UNC13C", "FCGBP", "SF3B1", "ELP1", "NES", "TRERF1")

# --- 1. Carica il file di dati ---
# Adatta questo percorso al tuo file reale
file_mutazioni_path <- "./resultsxmanuscript/Supplementary_table_2_BLCA_specific_gene_mutation_frequency_bmi_exon.csv" # <--- MODIFICA QUESTO CON IL NOME ESATTO DEL TUO FILE

# Prova a leggere come .xls/.xlsx. Se fallisce, prova come file di testo delimitato da tab.
dati_mutazioni <- tryCatch({
  read_excel(file_mutazioni_path)
}, error = function(e) {
  message("Errore nella lettura del file Excel, provando a leggere come file delimitato da tabulazioni...")
  read.delim(file_mutazioni_path, sep = "\t", header = TRUE) # Assumiamo tabulazioni come separatore
})

# --- 2. Preparazione dei Dati ---

# Rimuovi le righe con NA in bmi_category e filtra per i geni di interesse
dati_mutazioni_puliti <- dati_mutazioni %>%
  filter(!is.na(bmi_category)) %>%
  filter(gene %in% geni_di_interesse) %>%
  # Per la heatmap, vogliamo solo la presenza/assenza di mutazione per ogni paziente-gene
  # Usiamo distinct per contare ogni paziente una sola volta per gene, indipendentemente dal numero di mutazioni
  distinct(patients, gene, bmi_category, .keep_all = TRUE) %>% # Mantieni le altre colonne per riferimento se necessario
  mutate(has_mutation = 1) # Aggiungi una colonna che indica la presenza di mutazione

# Crea un dataframe completo con tutte le combinazioni gene-paziente
# Questo è necessario per mostrare anche le celle dove NON c'è mutazione (saranno bianche)
all_combinations <- expand_grid(
  patients = unique(dati_mutazioni_puliti$patients),
  gene = geni_di_interesse # Usa l'ordine definito
) %>%
  # Aggiungi bmi_category ai pazienti, usando un join
  left_join(dati_mutazioni_puliti %>% select(patients, bmi_category) %>% distinct(), by = "patients") %>%
  # Aggiungi le informazioni sulla mutazione
  left_join(dati_mutazioni_puliti %>% select(patients, gene, has_mutation), by = c("patients", "gene")) %>%
  mutate(has_mutation = ifelse(is.na(has_mutation), 0, has_mutation)) %>% # Riempi NA con 0
  filter(!is.na(bmi_category)) # Rimuovi i pazienti che non hanno una bmi_category valida

# Ordina i pazienti all'interno di ogni categoria BMI (sull'asse Y per la rotazione)
# Questo può essere fatto in modo casuale o basato su un altro criterio se disponibile
# Per ora, li ordiniamo alfabeticamente all'interno di ciascuna categoria per riproducibilità
all_combinations <- all_combinations %>%
  arrange(bmi_category, patients)

# Ordina i geni sull'asse X per la rotazione
# Mantieni l'ordine desiderato per i geni
all_combinations$gene <- factor(all_combinations$gene, levels = geni_di_interesse)

# Converti bmi_category in fattore con l'ordine desiderato per i facet E NUOVI NOMI
# L'ordine qui è CRITICO per il strip.background fill
all_combinations$bmi_category <- factor(all_combinations$bmi_category,
                                        levels = c("< 25", "≥ 25"), # Original values in the data
                                        labels = c("BMI < 25", "BMI ≥ 25")) # How they will appear on the plot


# --- 3. Creazione della Heatmap ---

p <- ggplot(all_combinations, aes(x = gene, y = patients, fill = as.factor(has_mutation))) + # SWAPPED X and Y for rotation
  geom_tile(color = "black", linewidth = 0.5) + # Aggiungi un bordo nero a ogni cella
  scale_fill_manual(
    values = c("0" = "white", "1" = "darkorchid4"), # 0 = nessuna mutazione (bianco), 1 = mutazione (viola scuro)
    labels = c("0" = "No Mutation", "1" = "Mutation Present"),
    name = NULL # Rimuovi il titolo della legenda
  ) +
  # Dividi per BMI usando facet_grid con bmi_category sulle colonne (vertical split)
  # scales = "free_y" e space = "free_y" permettono che le altezze dei pannelli siano diverse
  # se il numero di pazienti varia molto tra le categorie BMI.
  facet_grid(bmi_category ~ ., scales = "free_y", space = "free_y") + # BMI categories on Y-axis (vertical bands)
  labs(
#    title = "Heatmap of Nonsynonymous Mutations by Patient and Gene",
    x = "Gene", # Gene on X-axis (columns)
    y = "Patients" # Patients on Y-axis (rows)
  ) +
  theme_minimal() +
  theme(
    # Standardizzazione del Tema
    plot.title = element_text(hjust = 0.5, face = "bold", size = 64, family = "Arial", color = "black"),
    axis.title.x = element_text(size = 64, family = "Arial", color = "black", face = "bold", margin = margin(t = 20)),
    axis.title.y = element_text(size = 64, family = "Arial", color = "black", face = "bold", margin = margin(r = 20)),
    # Etichette dei geni sull'asse X, ruotate per leggibilità
    axis.text.x = element_text(angle = 45, hjust = 1, size = 69, family = "Arial", color = "black", face = "bold"),
    axis.ticks.x = element_blank(), # Keep tick marks visible with labels
    # Rimuovi le etichette dei singoli pazienti sull'asse Y per evitare sovrapposizioni
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), # Rimuovi i tick marks sull'asse Y
    legend.title = element_text(size = 57, family = "Arial", color = "black", face = "bold"),
    legend.text = element_text(size = 57, family = "Arial", color = "black", face = "bold"),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black"),
    strip.text = element_text(face = "bold", size = 57, family = "Arial", color = "black"),
    # *** FIX CRITICO QUI: Mappa esplicitamente le etichette ai colori ***
    strip.background = element_rect(fill = c("BMI < 25" = "orange3", "BMI ≥ 25" = "dodgerblue4"), color = "black")
  )

# --- 4. Salva il grafico in formato TIFF ---
output_file_heatmap <- "mutation_heatmap_bmi_stratified_rotated_elongated_final.tiff" # Nuovo nome file
ggsave(
  filename = output_file_heatmap,
  plot = p,
  width = 8, # Ridotto per mantenere proporzioni con l'aumento dell'altezza
  height = 12, # Aumentato per renderlo più "allungato" verticalmente
  dpi = 300,
  device = "tiff",
  bg = "white"
)

message(paste("Heatmap delle mutazioni generata e salvata come:", output_file_heatmap))