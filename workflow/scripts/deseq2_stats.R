

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
})

# ============================================================
# 1. Entrées et sorties Snakemake
# ============================================================

input_csv <- snakemake@input[["deseq2"]]

output_stat <- snakemake@output[["stat"]]

output_volcano_total <- snakemake@output[["volcano_total"]]
output_volcano_zoom  <- snakemake@output[["volcano_zoom"]]

output_go_total <- snakemake@output[["go_total"]]
output_go_up    <- snakemake@output[["go_enrich_up"]]
output_go_down  <- snakemake@output[["go_enrich_down"]]

output_deg_total <- snakemake@output[["deg_total"]]
output_deg_up    <- snakemake@output[["deg_up"]]
output_deg_down  <- snakemake@output[["deg_down"]]

output_directories <- unique(
  dirname(
    c(
      output_stat,
      output_volcano_total,
      output_volcano_zoom,
      output_deg_total,
      output_deg_up,
      output_deg_down,
      output_go_total,
      output_go_up,
      output_go_down
    )
  )
)

for (output_directory in output_directories) {
  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ============================================================
# 2. Lecture des résultats DESeq2
# ============================================================

resultats <- data.table::fread(
  input_csv,
  sep = ",",
  header = TRUE,
  data.table = FALSE,
  na.strings = c("NA", "NaN", "")
)

if (ncol(resultats) < 2) {
  stop(
    "Le fichier DESeq2 ne contient pas suffisamment de colonnes : ",
    input_csv
  )
}

# Le CSV produit par DESeq2 contient généralement les identifiants
# Ensembl dans la première colonne, parfois sans nom de colonne.
colnames(resultats)[1] <- "gene"

required_columns <- c(
  "gene",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)

missing_columns <- setdiff(
  required_columns,
  colnames(resultats)
)

if (length(missing_columns) > 0) {
  stop(
    "Colonnes manquantes dans le fichier DESeq2 : ",
    paste(missing_columns, collapse = ", ")
  )
}

# Conversion explicite des colonnes numériques.
numeric_columns <- c(
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)

resultats[numeric_columns] <- lapply(
  resultats[numeric_columns],
  as.numeric
)

# Nettoyage des identifiants.
resultats <- resultats %>%
  dplyr::mutate(
    gene = trimws(as.character(gene)),

    # ENSMUSG00000000001.5 -> ENSMUSG00000000001
    gene = sub("\\.[0-9]+$", "", gene)
  ) %>%
  dplyr::filter(
    !is.na(gene),
    gene != ""
  )

if (nrow(resultats) == 0) {
  stop(
    "Aucun identifiant de gène valide n'a été trouvé dans : ",
    input_csv
  )
}

# ============================================================
# 3. Annotation des identifiants Ensembl
# ============================================================

ensembl_keys <- unique(resultats$gene)

gene_info <- AnnotationDbi::select(
  x = org.Mm.eg.db,
  keys = ensembl_keys,
  keytype = "ENSEMBL",
  columns = c(
    "SYMBOL",
    "GENENAME",
    "ENTREZID"
  )
)

gene_info <- as.data.frame(
  gene_info,
  stringsAsFactors = FALSE
)

# AnnotationDbi::select() doit retourner la colonne utilisée comme clé.
# Vérification explicite afin d'éviter l'erreur :
# Error: object 'ENSEMBL' not found
if (!"ENSEMBL" %in% colnames(gene_info)) {
  stop(
    paste0(
      "La colonne ENSEMBL est absente du résultat de ",
      "AnnotationDbi::select(). Colonnes obtenues : ",
      paste(colnames(gene_info), collapse = ", ")
    )
  )
}

# Renommage en base R, sans évaluation non standard de dplyr.
colnames(gene_info)[
  colnames(gene_info) == "ENSEMBL"
] <- "gene"

colnames(gene_info)[
  colnames(gene_info) == "SYMBOL"
] <- "gene_symbol"

colnames(gene_info)[
  colnames(gene_info) == "GENENAME"
] <- "gene_name"

colnames(gene_info)[
  colnames(gene_info) == "ENTREZID"
] <- "entrez_id"

# Un identifiant Ensembl peut avoir plusieurs annotations.
# On conserve une seule ligne par identifiant pour éviter de dupliquer
# les lignes du tableau DESeq2 lors du left_join().
gene_info <- gene_info %>%
  dplyr::arrange(
    gene,
    is.na(gene_symbol),
    is.na(entrez_id),
    is.na(gene_name)
  ) %>%
  dplyr::distinct(
    gene,
    .keep_all = TRUE
  )

resultats <- resultats %>%
  dplyr::left_join(
    gene_info,
    by = "gene"
  ) %>%
  dplyr::mutate(
    gene_symbol = dplyr::if_else(
      is.na(gene_symbol) | gene_symbol == "",
      gene,
      gene_symbol
    )
  )

# ============================================================
# 4. Classification des gènes
# ============================================================

padj_threshold <- 0.05
log2fc_threshold <- 1

resultats <- resultats %>%
  dplyr::mutate(
    status = dplyr::case_when(
      !is.na(padj) &
        padj < padj_threshold &
        log2FoldChange >= log2fc_threshold ~ "Upregulated",

      !is.na(padj) &
        padj < padj_threshold &
        log2FoldChange <= -log2fc_threshold ~ "Downregulated",

      TRUE ~ "Not significant"
    ),

    # Valeur utilisée uniquement pour l'affichage du volcano plot.
    padj_plot = dplyr::case_when(
      is.na(padj) ~ 1,
      padj <= 0 ~ .Machine$double.xmin,
      TRUE ~ padj
    ),

    minus_log10_padj = -log10(padj_plot)
  )

# ============================================================
# 5. Tableau récapitulatif
# ============================================================

resultats <- resultats %>%
  dplyr::arrange(
    is.na(padj),
    padj,
    dplyr::desc(abs(log2FoldChange))
  )

data.table::fwrite(
  resultats,
  file = output_stat,
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

# ============================================================
# 6. Fonction de volcano plot
# ============================================================

create_volcano_plot <- function(
  df,
  title,
  x_limits = NULL,
  y_limits = NULL,
  max_labels = 20
) {

  label_df <- df %>%
    dplyr::filter(
      status != "Not significant",
      !is.na(padj),
      !is.na(log2FoldChange),
      is.finite(minus_log10_padj)
    ) %>%
    dplyr::arrange(
      padj,
      dplyr::desc(abs(log2FoldChange))
    ) %>%
    dplyr::slice_head(
      n = max_labels
    )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = log2FoldChange,
      y = minus_log10_padj,
      color = status
    )
  ) +
    ggplot2::geom_point(
      size = 1.5,
      alpha = 0.7,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Not significant" = "black",
        "Upregulated" = "red",
        "Downregulated" = "blue"
      ),
      breaks = c(
        "Upregulated",
        "Downregulated",
        "Not significant"
      )
    ) +
    ggplot2::geom_vline(
      xintercept = c(
        -log2fc_threshold,
        log2fc_threshold
      ),
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(padj_threshold),
      linetype = "dashed"
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(
        label = gene_symbol
      ),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.3,
      min.segment.length = 0,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        "BH-adjusted p-value < ",
        padj_threshold,
        " and |log2FC| ≥ ",
        log2fc_threshold
      ),
      x = "Log2 fold change",
      y = expression(-log[10]("BH-adjusted p-value")),
      color = "Status"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 16,
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10,
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(
        size = 13
      ),
      axis.text = ggplot2::element_text(
        size = 11
      ),
      legend.position = "right"
    )

  if (!is.null(x_limits) || !is.null(y_limits)) {
    p <- p +
      ggplot2::coord_cartesian(
        xlim = x_limits,
        ylim = y_limits
      )
  }

  return(p)
}

# ============================================================
# 7. Volcano plots
# ============================================================

comparison_name <- tools::file_path_sans_ext(
  basename(input_csv)
)

comparison_name <- sub(
  "_DESeq2_gene$",
  "",
  comparison_name
)

volcano_total <- create_volcano_plot(
  df = resultats,
  title = paste(
    "Volcano plot:",
    comparison_name
  ),
  max_labels = 20
)

volcano_zoom <- create_volcano_plot(
  df = resultats,
  title = paste(
    "Volcano plot zoom:",
    comparison_name
  ),
  x_limits = c(-5, 5),
  y_limits = c(0, 15),
  max_labels = 15
)

ggplot2::ggsave(
  filename = output_volcano_total,
  plot = volcano_total,
  width = 9,
  height = 7,
  dpi = 300
)

ggplot2::ggsave(
  filename = output_volcano_zoom,
  plot = volcano_zoom,
  width = 9,
  height = 7,
  dpi = 300
)

# ============================================================
# 8. Figure vide pour les enrichissements sans résultat
# ============================================================

save_empty_plot <- function(
  output_file,
  message_text
) {

  empty_plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      geom = "text",
      x = 0,
      y = 0,
      label = message_text,
      size = 5
    ) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_void()

  ggplot2::ggsave(
    filename = output_file,
    plot = empty_plot,
    width = 9,
    height = 7,
    dpi = 300
  )
}

# ============================================================
# 9. Univers des gènes pour GO
# ============================================================

# Tous les gènes testés par DESeq2.
gene_universe <- resultats %>%
  dplyr::filter(
    !is.na(gene),
    gene != ""
  ) %>%
  dplyr::pull(gene) %>%
  unique()

universe_annotation <- AnnotationDbi::select(
  x = org.Mm.eg.db,
  keys = gene_universe,
  keytype = "ENSEMBL",
  columns = "GO"
)

universe_annotation <- as.data.frame(
  universe_annotation,
  stringsAsFactors = FALSE
)

if (!"ENSEMBL" %in% colnames(universe_annotation)) {
  stop(
    paste0(
      "La colonne ENSEMBL est absente de l'annotation GO. ",
      "Colonnes obtenues : ",
      paste(colnames(universe_annotation), collapse = ", ")
    )
  )
}

if (!"GO" %in% colnames(universe_annotation)) {
  stop(
    paste0(
      "La colonne GO est absente de l'annotation GO. ",
      "Colonnes obtenues : ",
      paste(colnames(universe_annotation), collapse = ", ")
    )
  )
}

annotated_universe <- universe_annotation %>%
  dplyr::filter(
    !is.na(GO),
    GO != "",
    !is.na(ENSEMBL),
    ENSEMBL != ""
  ) %>%
  dplyr::pull(ENSEMBL) %>%
  unique()

# ============================================================
# 10. Fonction d'enrichissement GO
# ============================================================

run_go <- function(
  genes,
  output_file,
  title,
  show_category = 25
) {

  genes <- unique(
    genes[
      !is.na(genes) &
        genes != ""
    ]
  )

  # Seulement les gènes appartenant à l'univers annoté.
  genes <- intersect(
    genes,
    annotated_universe
  )

  if (length(genes) == 0) {
    save_empty_plot(
      output_file,
      paste(
        "Aucun gène admissible pour",
        title
      )
    )

    return(invisible(NULL))
  }

  ego <- clusterProfiler::enrichGO(
    gene = genes,
    universe = annotated_universe,
    OrgDb = org.Mm.eg.db,
    keyType = "ENSEMBL",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  ego_df <- as.data.frame(ego)

  if (nrow(ego_df) == 0) {
    save_empty_plot(
      output_file,
      paste(
        "Aucun terme GO significatif pour",
        title
      )
    )

    return(invisible(NULL))
  }

  go_plot <- enrichplot::dotplot(
    ego,
    showCategory = show_category
  ) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 15,
        hjust = 0.5
      ),
      axis.text.y = ggplot2::element_text(
        size = 9
      )
    )

  ggplot2::ggsave(
    filename = output_file,
    plot = go_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# ============================================================
# 11. Listes et tableaux de DEG
# ============================================================

deg_total_df <- resultats %>%
  dplyr::filter(
    !is.na(padj),
    padj < padj_threshold,
    !is.na(log2FoldChange),
    abs(log2FoldChange) >= log2fc_threshold
  ) %>%
  dplyr::arrange(
    padj,
    dplyr::desc(abs(log2FoldChange))
  )

deg_up_df <- resultats %>%
  dplyr::filter(
    !is.na(padj),
    padj < padj_threshold,
    !is.na(log2FoldChange),
    log2FoldChange >= log2fc_threshold
  ) %>%
  dplyr::arrange(
    padj,
    dplyr::desc(log2FoldChange)
  )

deg_down_df <- resultats %>%
  dplyr::filter(
    !is.na(padj),
    padj < padj_threshold,
    !is.na(log2FoldChange),
    log2FoldChange <= -log2fc_threshold
  ) %>%
  dplyr::arrange(
    padj,
    log2FoldChange
  )

# Vecteurs utilisés pour les enrichissements GO
deg_total <- unique(deg_total_df$gene)
deg_up <- unique(deg_up_df$gene)
deg_down <- unique(deg_down_df$gene)

# ============================================================
# 12. Export des DEG en TSV
# ============================================================

data.table::fwrite(
  deg_total_df,
  file = output_deg_total,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

data.table::fwrite(
  deg_up_df,
  file = output_deg_up,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

data.table::fwrite(
  deg_down_df,
  file = output_deg_down,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

# ============================================================
# 13. Enrichissements GO
# ============================================================

run_go(
  genes = deg_total,
  output_file = output_go_total,
  title = paste(
    comparison_name,
    "— all DEG"
  ),
  show_category = 30
)

run_go(
  genes = deg_up,
  output_file = output_go_up,
  title = paste(
    comparison_name,
    "— upregulated genes"
  ),
  show_category = 25
)

run_go(
  genes = deg_down,
  output_file = output_go_down,
  title = paste(
    comparison_name,
    "— downregulated genes"
  ),
  show_category = 25
)

# ============================================================
# 13. Résumé dans le log
# ============================================================

cat(
  "Comparison:", comparison_name, "\n",
  "Input file:", input_csv, "\n",
  "Genes in DESeq2 table:", nrow(resultats), "\n",
  "Genes mapped in org.Mm.eg.db:",
  sum(!is.na(resultats$entrez_id) | !is.na(resultats$gene_name)),
  "\n",
  "Significant DEG:", length(deg_total), "\n",
  "Upregulated:", length(deg_up), "\n",
  "Downregulated:", length(deg_down), "\n",
  "GO universe:", length(annotated_universe), "\n"
)
