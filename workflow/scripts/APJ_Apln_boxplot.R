#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(openxlsx)
  library(ragg)
  library(svglite)
})

# ============================================================
# 1. Paramètres
# ============================================================

dds_file <- "results/deseq2_WT_only/dds_WT_only.rds"
out_dir  <- "results/boxplots"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

font_family <- "Liberation Sans"

genes_of_interest <- data.frame(
  gene_id   = c("ENSMUSG00000044338", "ENSMUSG00000037010"),
  gene_name = c("Aplnr", "Apln"),
  stringsAsFactors = FALSE
)

# ============================================================
# 2. Chargement du DDS WT-only
# ============================================================

dds <- readRDS(dds_file)

cat("\nSamples présents dans le DDS :\n")
print(colnames(dds))

cat("\nConditions :\n")
print(table(dds$condition))

if (ncol(dds) != 8) {
  warning(
    "Le DDS contient ",
    ncol(dds),
    " échantillons au lieu des 8 WT attendus."
  )
}

# ============================================================
# 3. Résultat DESeq2 global
# ============================================================

res <- results(
  dds,
  contrast = c("condition", "patho", "chow")
)

# ============================================================
# 4. Fonction significativité
# ============================================================

get_sig_label <- function(padj) {
  if (is.na(padj)) {
    return("ns")
  } else if (padj < 0.0001) {
    return("****")
  } else if (padj < 0.001) {
    return("***")
  } else if (padj < 0.01) {
    return("**")
  } else if (padj < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# ============================================================
# 5. Extraction des données pour chaque gène
# ============================================================

plot_list  <- list()
stats_list <- list()

for (i in seq_len(nrow(genes_of_interest))) {

  gene_id   <- genes_of_interest$gene_id[i]
  gene_name <- genes_of_interest$gene_name[i]

  cat("\n============================================\n")
  cat("Traitement de :", gene_name, "(", gene_id, ")\n")
  cat("============================================\n")

  if (!gene_id %in% rownames(res)) {
    warning("Le gène ", gene_id, " n'est pas présent dans le DDS. Il sera ignoré.")
    next
  }

  gene_res <- res[gene_id, ]
  padj <- gene_res$padj
  sig_label <- get_sig_label(padj)

  print(gene_res)

  # normalized counts
  df <- plotCounts(
    dds,
    gene = gene_id,
    intgroup = "condition",
    returnData = TRUE
  )

  df$sample <- rownames(df)

  df$condition <- factor(
    df$condition,
    levels = c("chow", "patho"),
    labels = c("Chow", "HFpEF")
  )

  # positions rapprochées
  df$x_pos <- ifelse(df$condition == "Chow", 1.20, 1.80)

  # transformation
  df$log2_count <- log2(df$count + 1)

  # annotations gène
  df$gene_id   <- gene_id
  df$gene_name <- gene_name
  df$padj      <- padj
  df$sig_label <- sig_label

  plot_list[[gene_name]] <- df

  stats_list[[gene_name]] <- data.frame(
    gene = gene_name,
    gene_id = gene_id,
    comparison = "HFpEF WT vs Chow WT",
    baseMean = gene_res$baseMean,
    log2FoldChange = gene_res$log2FoldChange,
    lfcSE = gene_res$lfcSE,
    stat = gene_res$stat,
    pvalue = gene_res$pvalue,
    padj = gene_res$padj,
    sig_label = sig_label,
    row.names = NULL
  )

  cat(
    "FDR = ",
    signif(padj, 4),
    " -> ",
    sig_label,
    "\n",
    sep = ""
  )
}

plot_df  <- do.call(rbind, plot_list)
stats_df <- do.call(rbind, stats_list)

# ============================================================
# 6. Tableau export
# ============================================================

df_export <- data.frame(
  gene = plot_df$gene_name,
  gene_id = plot_df$gene_id,
  sample = plot_df$sample,
  count = plot_df$count,
  condition = as.character(plot_df$condition),
  log2_count = plot_df$log2_count,
  row.names = NULL
)

cat("\nValeurs utilisées pour la figure :\n")
print(df_export)

# ============================================================
# 7. Export XLSX
# ============================================================

xlsx_file <- file.path(
  out_dir,
  "Aplnr_Apln_WT_normalized_counts.xlsx"
)

wb <- createWorkbook()

addWorksheet(wb, "Normalized_counts")
writeData(
  wb,
  sheet = "Normalized_counts",
  x = df_export
)

setColWidths(
  wb,
  sheet = "Normalized_counts",
  cols = 1:ncol(df_export),
  widths = "auto"
)

addWorksheet(wb, "DESeq2_statistics")
writeData(
  wb,
  sheet = "DESeq2_statistics",
  x = stats_df
)

setColWidths(
  wb,
  sheet = "DESeq2_statistics",
  cols = 1:ncol(stats_df),
  widths = "auto"
)

saveWorkbook(
  wb,
  xlsx_file,
  overwrite = TRUE
)

# ============================================================
# 8. Position des 4 groupes sur un seul graphique
# ============================================================

plot_df$gene_name <- factor(
  plot_df$gene_name,
  levels = c("Aplnr", "Apln")
)

plot_df$group <- interaction(
  plot_df$gene_name,
  plot_df$condition,
  sep = "_"
)

# Positions :
# Aplnr : 1.00 / 1.45
# Apln  : 2.15 / 2.60

plot_df$x_pos <- NA_real_

plot_df$x_pos[
  plot_df$gene_name == "Aplnr" &
  plot_df$condition == "Chow"
] <- 1.00

plot_df$x_pos[
  plot_df$gene_name == "Aplnr" &
  plot_df$condition == "HFpEF"
] <- 1.45

plot_df$x_pos[
  plot_df$gene_name == "Apln" &
  plot_df$condition == "Chow"
] <- 2.15

plot_df$x_pos[
  plot_df$gene_name == "Apln" &
  plot_df$condition == "HFpEF"
] <- 2.60


# ============================================================
# 9. Significativité + paramètres axe Y
# ============================================================

apl_nr_sig <- stats_df$sig_label[
  stats_df$gene == "Aplnr"
]

apln_sig <- stats_df$sig_label[
  stats_df$gene == "Apln"
]

# Axe Y
plot_ymin <- 8.80
plot_ymax <- 11.80

y_breaks <- seq(
  from = 8.80,
  to = 11.80,
  by = 0.50
)

# Hauteur automatique des crochets
y_max <- max(plot_df$log2_count, na.rm = TRUE)
y_min <- min(plot_df$log2_count, na.rm = TRUE)

y_range <- y_max - y_min

if (y_range == 0) {
  y_range <- 1
}

bracket_Aplnr <- y_max + 0.08 * y_range
bracket_Apln  <- y_max + 0.08 * y_range

tip_length <- 0.04 * y_range
text_y <- y_max + 0.13 * y_range


# ============================================================
# 10. Figure
# ============================================================

p <- ggplot(
  plot_df,
  aes(
    x = x_pos,
    y = log2_count,
    fill = condition
  )
) +

  geom_boxplot(
    aes(
      group = group,
      colour = condition
    ),
    width = 0.22,
    outlier.shape = NA,
    linewidth = 0.9
  ) +

  geom_jitter(
    aes(
      colour = condition,
      group = group
    ),
    width = 0.025,
    height = 0,
    size = 3.0,
    alpha = 1
  ) +

  annotate(
    "segment",
    x = 1.00,
    xend = 1.45,
    y = bracket_Aplnr,
    yend = bracket_Aplnr,
    linewidth = 0.8
  ) +

  annotate(
    "segment",
    x = 1.00,
    xend = 1.00,
    y = bracket_Aplnr,
    yend = bracket_Aplnr - tip_length,
    linewidth = 0.8
  ) +

  annotate(
    "segment",
    x = 1.45,
    xend = 1.45,
    y = bracket_Aplnr,
    yend = bracket_Aplnr - tip_length,
    linewidth = 0.8
  ) +

  annotate(
    "text",
    x = 1.225,
    y = text_y,
    label = apl_nr_sig,
    family = font_family,
    fontface = "bold",
    size = 5
  ) +

  annotate(
    "segment",
    x = 2.15,
    xend = 2.60,
    y = bracket_Apln,
    yend = bracket_Apln,
    linewidth = 0.8
  ) +

  annotate(
    "segment",
    x = 2.15,
    xend = 2.15,
    y = bracket_Apln,
    yend = bracket_Apln - tip_length,
    linewidth = 0.8
  ) +

  annotate(
    "segment",
    x = 2.60,
    xend = 2.60,
    y = bracket_Apln,
    yend = bracket_Apln - tip_length,
    linewidth = 0.8
  ) +

  annotate(
    "text",
    x = 2.375,
    y = text_y,
    label = apln_sig,
    family = font_family,
    fontface = "bold",
    size = 5
  ) +

  scale_fill_manual(
    values = c(
      "Chow"  = "#A0A0A440",
      "HFpEF" = "#FF808040"
    )
  ) +

  scale_colour_manual(
    values = c(
      "Chow"  = "#A0A0A4",
      "HFpEF" = "#FF8080"
    )
  ) +

  # ==========================================================
  # Axe X : 4 boxplots
  # ==========================================================

  scale_x_continuous(
    limits = c(0.70, 2.90),

    breaks = c(
      1.00,
      1.45,
      2.15,
      2.60
    ),

    labels = c(
      "Chow",
      "HFpEF",
      "Chow",
      "HFpEF"
    ),

    expand = c(0, 0)
  ) +

  # ==========================================================
  # Axe Y
  # ==========================================================

  scale_y_continuous(
    breaks = y_breaks,
    expand = expansion(mult = c(0, 0))
  ) +

  labs(
    x = NULL,
    y = expression(
      log[2] * "(DESeq2 normalized counts + 1)"
    )
  ) +

  coord_cartesian(
    ylim = c(plot_ymin, plot_ymax),
    clip = "off"
  ) +

  # ==========================================================
  # Style
  # ==========================================================

  theme_classic(
    base_size = 14,
    base_family = font_family
  ) +

  theme(
    legend.position = "none",

    axis.title.y = element_text(
      face = "bold",
      size = 13,
      colour = "black"
    ),

    axis.text.x = element_text(
      face = "bold",
      size = 12,
      colour = "black"
    ),

    axis.text.y = element_text(
      face = "bold",
      size = 11,
      colour = "black"
    ),

    axis.line = element_line(
      linewidth = 0.9,
      colour = "black"
    ),

    axis.ticks = element_line(
      linewidth = 0.9,
      colour = "black"
    ),

    plot.margin = margin(
      t = 20,
      r = 10,
      b = 35,
      l = 10
    )
  ) +

  # ==========================================================
  # Nom des deux gènes sous les groupes
  # ==========================================================

  annotate(
    "text",
    x = 1.225,
    y = -Inf,
    label = "Aplnr",
    vjust = 4.0,
    fontface = "bold",
    family = font_family,
    size = 4.5
  ) +

  annotate(
    "text",
    x = 2.375,
    y = -Inf,
    label = "Apln",
    vjust = 4.0,
    fontface = "bold",
    family = font_family,
    size = 4.5
  )
  
# ============================================================
# 10. Fichiers de sortie
# ============================================================

png_file <- file.path(
  out_dir,
  "Aplnr_Apln_WT_boxplot.png"
)

pdf_file <- file.path(
  out_dir,
  "Aplnr_Apln_WT_boxplot.pdf"
)

svg_file <- file.path(
  out_dir,
  "Aplnr_Apln_WT_boxplot.svg"
)

# ============================================================
# 11. Export PNG
# ============================================================

ggsave(
  filename = png_file,
  plot = p,
  device = ragg::agg_png,
  width = 6.4,
  height = 4.8,
  units = "in",
  res = 600,
  bg = "white"
)

# ============================================================
# 12. Export PDF
# ============================================================

ggsave(
  filename = pdf_file,
  plot = p,
  device = grDevices::cairo_pdf,
  width = 6.4,
  height = 4.8,
  units = "in",
  bg = "white"
)

# ============================================================
# 13. Export SVG
# ============================================================

ggsave(
  filename = svg_file,
  plot = p,
  device = svglite::svglite,
  width = 6.4,
  height = 4.8,
  units = "in",
  bg = "white"
)

# ============================================================
# 14. Fin
# ============================================================

cat(
  "\n============================================\n",
  "Analyse terminée\n",
  "============================================\n\n",
  "PNG  : ", png_file, "\n",
  "PDF  : ", pdf_file, "\n",
  "SVG  : ", svg_file, "\n",
  "XLSX : ", xlsx_file, "\n",
  sep = ""
)
