library(ggplot2)
library(ggrepel)
library(dplyr)
library(showtext)

# Attiva il rendering dei font
showtext_auto()

# Lista dei geni da etichettare
cancer_geni_da_etichettare <- c("blca_BRCA2", "blca_DNAH9", "blca_GRIA4", "blca_PLXNA4", "blca_UNC13C",
                                "blca_FCGBP", "blca_SF3B1", "blca_ELP1", "blca_NES", "blca_TRERF1")

# Funzione aggiornata con tutte le correzioni
crea_volcano_plot_da_dati_corretto <- function(df, nome_file_output = "volcano_plot_da_file.tiff",
                                               x_col = "estimate", padj_col = "p.adj", cancer_gene_col = "cancer_gene",
                                               p_val_threshold = 0.05,
                                               title = "Volcano Plot da File",
                                               cancer_geni_etichetta = cancer_geni_da_etichettare) {
  
  # Controllo iniziale sul dataframe
  if (is.null(df) || nrow(df) == 0) {
    message("Il dataframe fornito è vuoto. Nessun volcano plot generato.")
    return(NULL)
  }
  
  # Converti le colonne necessarie a numeriche
  if (!all(c(x_col, padj_col, cancer_gene_col) %in% names(df))) {
    stop("Le colonne 'estimate', 'p.adj' o 'cancer_gene' non sono presenti nel dataframe.")
  }
  
  df[[x_col]] <- as.numeric(df[[x_col]])
  df[[padj_col]] <- as.numeric(df[[padj_col]])
  
  # Identifica i geni significativi
  df <- df %>%
    mutate(significant = ifelse(.data[[padj_col]] < p_val_threshold, "yes", "no"))
  
  # Crea colonna per le etichette
  df <- df %>%
    mutate(etichetta = ifelse(.data[[cancer_gene_col]] %in% cancer_geni_etichetta,
                              .data[[cancer_gene_col]], NA_character_))
  
  # Calcola i limiti dell'asse x con più spazio
  x_min <- min(df[[x_col]], na.rm = TRUE)
  x_max <- max(df[[x_col]], na.rm = TRUE)
  x_range <- max(abs(x_min), abs(x_max))
  x_limits <- c(-x_range * 1.35, x_range * 1.35)
  
  # Prepara sottotitolo (commentato, ma mantenuto per riferimento)
  # num_sig <- sum(df$significant == "yes", na.rm = TRUE)
  # plot_subtitle <- paste(num_sig, "significant genes (p.adj <", p_val_threshold, ")")
  
  # Crea il plot
  p <- ggplot(df, aes(x = !!sym(x_col), y = -log10(!!sym(padj_col)), color = significant)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_vline(xintercept = 0, color = "black", linetype = "solid") +
    geom_hline(yintercept = -log10(p_val_threshold), color = "blue", linetype = "dashed", alpha = 0.7) +
    
    # Unico geom_text_repel per tutte le etichette
    ggrepel::geom_text_repel(
      aes(label = etichetta),
      size = 28,
      color = "black",
      na.rm = TRUE,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.05,
      segment.size = 0.8,
      min.segment.length = 0,
      force = 4,
      fontface = "bold",
      nudge_x = 0.1,
      nudge_y = -0.2,
      segment.color = "grey50",
      segment.alpha = 0.6,
      family = "Arial",
      seed = 42
    ) +
    
    # *** Modifica qui per il colore viola (darkorchid4) per i punti significativi ***
    scale_color_manual(
      values = c("no" = "grey70", "yes" = "darkorchid4"), # "darkorchid4" for significant points
      labels = c("no" = "Non significant", "yes" = "Significant"),
      name = NULL
    ) +
    
    labs(
      x = "Effect size (estimate)",
      y = "-log10(Adjusted P-value)"
      # title = title,
      # subtitle = plot_subtitle
    ) +
    
    xlim(x_limits) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      # plot.title = element_text(hjust = 0.5, face = "bold", size = 64, family = "Arial", color = "black"),
      # plot.subtitle = element_text(hjust = 0.5, size = 50, family = "Arial", color = "black", margin = margin(b = 20)),
      axis.title = element_text(size = 70, family = "Arial", color = "black", face = "bold"),
      axis.text = element_text(size = 75, family = "Arial", color = "black", face = "bold"),
      legend.text = element_text(size = 65, family = "Arial", color = "black", face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black"),
      plot.margin = margin(1, 2, 1, 1, "cm")
    )
  
  # Salva il plot
  ggsave(nome_file_output, plot = p, width = 12, height = 8, units = "in", dpi = 300, bg = "white")
  
  message(paste("Volcano plot generato e salvato come:", nome_file_output))
  return(p)
}

# --- Caricamento dati e generazione plot ---
file_path <- "volcano_ns_bmi_data_filtered_sig.csv" # Assicurati che questo percorso sia corretto
dati_caricati <- tryCatch({
  read.csv(file_path, sep = ",")
}, error = function(e) {
  stop(paste("Errore nella lettura del file:", file_path, "-", e$message))
})

if (!is.null(dati_caricati)) {
  volcano_plot <- crea_volcano_plot_da_dati_corretto(
    df = dati_caricati,
    nome_file_output = "volcano_ns_bmi_final_segments_touching_purple.tiff", # Nuovo nome file per chiarezza
    title = "Volcano Plot - BMI Association with Mutation Presence"
  )
}