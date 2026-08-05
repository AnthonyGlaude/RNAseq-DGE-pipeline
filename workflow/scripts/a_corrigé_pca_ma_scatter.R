

# ============================================================
# 0. Packages
# ============================================================
#
# module load r/4.6.1
# Rscript ce_script.R
#
# ============================================================

required_packages <- c(
  "DESeq2",
  "ggplot2",
  "ggrepel",
  "plotly",
  "htmlwidgets"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {

  stop(
    paste0(
      "Packages R manquants : ",
      paste(
        missing_packages,
        collapse = ", "
      ),
      "\n\nInstalle-les avant de lancer le script."
    )
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(plotly)
  library(htmlwidgets)
})

# ============================================================
# 1. Fonction de sauvegarde des widgets HTML
# ============================================================
#
# Le script essaie d'abord de produire un seul fichier HTML
# autonome avec selfcontained = TRUE.
#
# Si Pandoc n'est pas disponible, il utilise automatiquement
# selfcontained = FALSE. Dans ce cas, un dossier contenant les
# dépendances JavaScript sera créé à côté du fichier HTML.
#
# ============================================================

save_interactive_widget <- function(
  widget,
  file
) {

  try_selfcontained <- tryCatch(

    {

      htmlwidgets::saveWidget(
        widget = widget,
        file = file,
        selfcontained = TRUE
      )

      TRUE
    },

    error = function(e) {

      message(
        "Sauvegarde HTML autonome impossible : ",
        conditionMessage(e)
      )

      FALSE
    }
  )

  if (!try_selfcontained) {

    message(
      "Nouvelle tentative avec selfcontained = FALSE."
    )

    htmlwidgets::saveWidget(
      widget = widget,
      file = file,
      selfcontained = FALSE
    )
  }

  message(
    "Figure interactive sauvegardée : ",
    file
  )
}

# ============================================================
# 2. Chemins et paramètres
# ============================================================

# 
script_dir <- "/home/glaudea/scratch/glaudea/test_souris_rnaseq/RNA_seq-analysis/workflow/results/deseq2"

dds_file <- file.path(
  script_dir,
  "dds.rds"
)

pca_file <- file.path(
  script_dir,
  "PCA_global.png"
)

ma_dir <- file.path(
  script_dir,
  "MAplots"
)

scatter_dir <- file.path(
  script_dir,
  "Scatterplots"
)

interactive_dir <- file.path(
  script_dir,
  "Interactive"
)

interactive_ma_dir <- file.path(
  interactive_dir,
  "MAplots"
)

interactive_scatter_dir <- file.path(
  interactive_dir,
  "Scatterplots"
)

dir.create(
  ma_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  scatter_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  interactive_ma_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  interactive_scatter_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Nombre de gènes utilisés pour la PCA
ntop <- 500

# Seuils pour identifier les DEG
padj_threshold <- 0.05
lfc_threshold <- 1

# Limites verticales affichées pour les MA plots
ma_ylim <- c(
  -5,
  5
)

message("DDS utilisé : ", dds_file)
message("PCA : ", pca_file)
message("MA plots PNG : ", ma_dir)
message("Scatterplots PNG : ", scatter_dir)
message("MA plots interactifs : ", interactive_ma_dir)
message(
  "Scatterplots interactifs : ",
  interactive_scatter_dir
)

# ============================================================
# 3. Lecture du DESeqDataSet
# ============================================================

if (!file.exists(dds_file)) {

  stop(
    "Fichier introuvable : ",
    dds_file
  )
}

dds <- readRDS(
  dds_file
)

metadata <- as.data.frame(
  colData(dds)
)

message("Métadonnées disponibles :")
print(metadata)

message("Colonnes disponibles :")
print(
  colnames(metadata)
)

# ============================================================
# 4. Création du facteur group
# ============================================================
#
# Le facteur group doit contenir exactement :
#
# Control_WT
# Control_KO
# HFpEF_WT
# HFpEF_KO
#
# Le script cherche d'abord une colonne contenant directement
# ces quatre groupes.
#
# Sinon, il construit group avec :
#
# condition + genotype
#
# ============================================================

expected_groups <- c(
  "Control_WT",
  "Control_KO",
  "HFpEF_WT",
  "HFpEF_KO"
)

group_column <- NULL

for (column_name in colnames(metadata)) {

  column_values <- unique(
    as.character(
      metadata[[column_name]]
    )
  )

  if (
    all(
      expected_groups %in% column_values
    )
  ) {

    group_column <- column_name
    break
  }
}

if (!is.null(group_column)) {

  message(
    "Les quatre groupes ont été trouvés dans la colonne : ",
    group_column
  )

  colData(dds)$group <- factor(
    as.character(
      colData(dds)[[group_column]]
    ),
    levels = expected_groups
  )

} else if (
  all(
    c(
      "condition",
      "genotype"
    ) %in% colnames(metadata)
  )
) {

  message(
    paste0(
      "Création de group à partir des colonnes ",
      "'condition' et 'genotype'."
    )
  )

  condition_clean <- trimws(
    as.character(
      colData(dds)$condition
    )
  )

  genotype_clean <- trimws(
    as.character(
      colData(dds)$genotype
    )
  )

  # ----------------------------------------------------------
  # Uniformisation de la condition
  # ----------------------------------------------------------

  condition_clean[
    tolower(condition_clean) %in%
      c(
        "chow",
        "control",
        "controle",
        "contrôle",
        "normal"
      )
  ] <- "Control"

  condition_clean[
    grepl(
      pattern = "hfpef|patho",
      x = condition_clean,
      ignore.case = TRUE
    )
  ] <- "HFpEF"

  # ----------------------------------------------------------
  # Uniformisation du génotype
  # ----------------------------------------------------------

  genotype_clean[
    grepl(
      pattern = "^wt$|wild",
      x = genotype_clean,
      ignore.case = TRUE
    )
  ] <- "WT"

  genotype_clean[
    grepl(
      pattern = "^ko$|knockout|knock-out",
      x = genotype_clean,
      ignore.case = TRUE
    )
  ] <- "KO"

  colData(dds)$group <- factor(
    paste(
      condition_clean,
      genotype_clean,
      sep = "_"
    ),
    levels = expected_groups
  )

} else {

  stop(
    paste0(
      "Impossible de créer les quatre groupes.\n\n",
      "Le tableau de métadonnées doit contenir soit une colonne ",
      "avec les valeurs :\n",
      paste(
        expected_groups,
        collapse = ", "
      ),
      "\n\nou deux colonnes nommées condition et genotype."
    )
  )
}

message("Répartition finale des groupes :")

print(
  table(
    colData(dds)$group,
    useNA = "ifany"
  )
)

if (
  any(
    is.na(
      colData(dds)$group
    )
  )
) {

  metadata_problem <- as.data.frame(
    colData(dds)
  )

  message(
    "Échantillons qui n'ont pas pu être assignés :"
  )

  print(
    metadata_problem[
      is.na(metadata_problem$group),
      ,
      drop = FALSE
    ]
  )

  stop(
    paste0(
      "Au moins un échantillon n'a pas pu être assigné à ",
      "Control_WT, Control_KO, HFpEF_WT ou HFpEF_KO."
    )
  )
}

missing_groups <- setdiff(
  expected_groups,
  unique(
    as.character(
      colData(dds)$group
    )
  )
)

if (length(missing_groups) > 0) {

  stop(
    "Groupes absents du DDS : ",
    paste(
      missing_groups,
      collapse = ", "
    )
  )
}

# ============================================================
# 5.  DESeq2
# ============================================================
#
# 
#
# ~ group
#
# Chaque contraste comparera ensuite deux niveaux de group.
#
# ============================================================

design(dds) <- ~ group

dds <- DESeq(
  dds,
  quiet = FALSE
)

message("Coefficients disponibles :")

print(
  resultsNames(dds)
)

# 
corrected_dds_file <- file.path(
  script_dir,
  "dds_group_design.rds"
)

saveRDS(
  dds,
  corrected_dds_file
)

message(
  "DDS avec design ~ group sauvegardé : ",
  corrected_dds_file
)

# ============================================================
# 6. Transformation VST pour la PCA
# ============================================================

vsd <- vst(
  dds,
  blind = FALSE
)

vst_matrix <- assay(
  vsd
)

gene_variances <- apply(
  vst_matrix,
  1,
  var
)

selected_genes <- order(
  gene_variances,
  decreasing = TRUE
)[
  seq_len(
    min(
      ntop,
      length(gene_variances)
    )
  )
]

vst_top <- vst_matrix[
  selected_genes,
  ,
  drop = FALSE
]

# ============================================================
# 7. PCA
# ============================================================

pca <- prcomp(
  t(vst_top),
  center = TRUE,
  scale. = FALSE
)

pca_df <- as.data.frame(
  pca$x
)

pca_df$sample <- rownames(
  pca_df
)

metadata <- as.data.frame(
  colData(dds)
)

metadata$sample <- rownames(
  metadata
)

# match() conserve exactement l'ordre des échantillons

metadata <- metadata[
  match(
    pca_df$sample,
    metadata$sample
  ),
  ,
  drop = FALSE
]

pca_df <- cbind(
  pca_df,
  metadata[
    ,
    setdiff(
      colnames(metadata),
      "sample"
    ),
    drop = FALSE
  ]
)

percent_var <- round(
  100 * pca$sdev^2 /
    sum(pca$sdev^2),
  digits = 1
)

p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = group,
    label = sample
  )
) +
  geom_point(
    size = 4,
    alpha = 0.85
  ) +
  geom_text_repel(
    size = 3.5,
    show.legend = FALSE
  ) +
  labs(
    title = paste0(
      "PCA : ",
      ntop,
      " gènes les plus variables"
    ),
    x = paste0(
      "PC1 (",
      percent_var[1],
      " %)"
    ),
    y = paste0(
      "PC2 (",
      percent_var[2],
      " %)"
    ),
    color = "Groupe"
  ) +
  theme_minimal(
    base_size = 14
  ) +
  theme(
    legend.position = "top",
    plot.title = element_text(
      hjust = 0.5
    )
  )

ggsave(
  filename = pca_file,
  plot = p_pca,
  width = 9,
  height = 7,
  dpi = 300
)

message(
  "PCA sauvegardée : ",
  pca_file
)

# ============================================================
# 8. Comparaisons
# ============================================================
#
# Le premier groupe est le numérateur.
#
# Exemple :
#
# Control_KO / Control_WT
#
# log2FoldChange positif :
# expression plus élevée dans Control_KO.
#
# log2FoldChange négatif :
# expression plus élevée dans Control_WT.
#
# ============================================================

comparisons <- list(

  Control_KO_vs_Control_WT = c(
    "Control_KO",
    "Control_WT"
  ),

  HFpEF_KO_vs_HFpEF_WT = c(
    "HFpEF_KO",
    "HFpEF_WT"
  ),

  HFpEF_WT_vs_Control_WT = c(
    "HFpEF_WT",
    "Control_WT"
  ),

  HFpEF_KO_vs_Control_KO = c(
    "HFpEF_KO",
    "Control_KO"
  )
)

# ============================================================
# 9. Comptes normalisés pour les scatterplots
# ============================================================

normalized_counts <- counts(
  dds,
  normalized = TRUE
)

sample_groups <- as.character(
  colData(dds)$group
)

# ============================================================
# 10. MA plots et scatterplots
# ============================================================

summary_list <- list()

for (
  comparison_name in names(comparisons)
) {

  groups <- comparisons[[comparison_name]]

  numerator_group <- groups[1]
  denominator_group <- groups[2]

  message(
    "=================================================="
  )

  message(
    "Analyse : ",
    numerator_group,
    " vs ",
    denominator_group
  )

  # ----------------------------------------------------------
  # Résultats DESeq2
  # ----------------------------------------------------------

  res <- results(
    dds,
    contrast = c(
      "group",
      numerator_group,
      denominator_group
    ),
    alpha = padj_threshold
  )

  res_df <- as.data.frame(
    res
  )

  res_df$gene <- rownames(
    res_df
  )

  # ----------------------------------------------------------
  # Catégorie différentielle
  # ----------------------------------------------------------

  res_df$status <- "Non significatif"

  res_df$status[
    !is.na(res_df$padj) &
      res_df$padj < padj_threshold &
      !is.na(res_df$log2FoldChange) &
      res_df$log2FoldChange >= lfc_threshold
  ] <- paste0(
    "Plus élevé dans ",
    numerator_group
  )

  res_df$status[
    !is.na(res_df$padj) &
      res_df$padj < padj_threshold &
      !is.na(res_df$log2FoldChange) &
      res_df$log2FoldChange <= -lfc_threshold
  ] <- paste0(
    "Plus élevé dans ",
    denominator_group
  )

  status_levels <- c(
    "Non significatif",
    paste0(
      "Plus élevé dans ",
      numerator_group
    ),
    paste0(
      "Plus élevé dans ",
      denominator_group
    )
  )

  res_df$status <- factor(
    res_df$status,
    levels = status_levels
  )

  # ----------------------------------------------------------
  # Texte affiché lors du survol du MA plot
  # ----------------------------------------------------------

  res_df$hover_text <- paste0(
    "<b>Gène : ",
    res_df$gene,
    "</b>",
    "<br>BaseMean : ",
    format(
      round(
        res_df$baseMean,
        3
      ),
      big.mark = " ",
      scientific = FALSE
    ),
    "<br>log2FC : ",
    round(
      res_df$log2FoldChange,
      3
    ),
    "<br>p-value : ",
    format(
      res_df$pvalue,
      scientific = TRUE,
      digits = 3
    ),
    "<br>padj : ",
    format(
      res_df$padj,
      scientific = TRUE,
      digits = 3
    ),
    "<br>Statut : ",
    res_df$status
  )

  # ----------------------------------------------------------
  # Sauvegarde du tableau DESeq2
  # ----------------------------------------------------------

  results_file <- file.path(
    script_dir,
    paste0(
      "DESeq2_",
      comparison_name,
      ".csv"
    )
  )

  write.csv(
    res_df,
    results_file,
    row.names = FALSE
  )

  message(
    "Résultats DESeq2 sauvegardés : ",
    results_file
  )

  # ==========================================================
  # 10.1 MA plot statique
  # ==========================================================

  ma_file <- file.path(
    ma_dir,
    paste0(
      "MAplot_",
      comparison_name,
      ".png"
    )
  )

  png(
    filename = ma_file,
    width = 1800,
    height = 1400,
    res = 200
  )

  plotMA(
    res,
    alpha = padj_threshold,
    ylim = ma_ylim,
    main = gsub(
      "_",
      " ",
      comparison_name
    )
  )

  abline(
    h = c(
      -lfc_threshold,
      lfc_threshold
    ),
    lty = 2
  )

  dev.off()

  message(
    "MA plot statique sauvegardé : ",
    ma_file
  )

  # ==========================================================
  # 10.2 MA plot interactif
  # ==========================================================

  ma_interactive_df <- res_df[
    is.finite(res_df$baseMean) &
      is.finite(res_df$log2FoldChange),
    ,
    drop = FALSE
  ]

  # log10(baseMean + 1) évite les problèmes avec baseMean = 0

  ma_interactive_df$log10_baseMean <- log10(
    ma_interactive_df$baseMean + 1
  )

  p_ma_interactive <- plot_ly(
    data = ma_interactive_df,
    x = ~log10_baseMean,
    y = ~log2FoldChange,
    color = ~status,
    colors = "Set1",
    type = "scatter",
    mode = "markers",
    text = ~hover_text,
    hoverinfo = "text",
    marker = list(
      size = 5,
      opacity = 0.55
    )
  )

  p_ma_interactive <- p_ma_interactive |>
    layout(
      title = list(
        text = gsub(
          "_",
          " ",
          comparison_name
        ),
        x = 0.5
      ),
      xaxis = list(
        title = "log10(baseMean + 1)",
        zeroline = FALSE
      ),
      yaxis = list(
        title = "log2 fold change",
        range = ma_ylim,
        zeroline = TRUE
      ),
      legend = list(
        title = list(
          text = "Résultat DESeq2"
        ),
        orientation = "h",
        x = 0,
        y = 1.12
      ),
      shapes = list(

        # Ligne centrale log2FC = 0
        list(
          type = "line",
          x0 = 0,
          x1 = 1,
          xref = "paper",
          y0 = 0,
          y1 = 0,
          line = list(
            dash = "solid",
            width = 1
          )
        ),

        # Seuil log2FC positif
        list(
          type = "line",
          x0 = 0,
          x1 = 1,
          xref = "paper",
          y0 = lfc_threshold,
          y1 = lfc_threshold,
          line = list(
            dash = "dash",
            width = 1
          )
        ),

        # Seuil log2FC négatif
        list(
          type = "line",
          x0 = 0,
          x1 = 1,
          xref = "paper",
          y0 = -lfc_threshold,
          y1 = -lfc_threshold,
          line = list(
            dash = "dash",
            width = 1
          )
        )
      ),
      hovermode = "closest"
    ) |>
    config(
      displaylogo = FALSE,
      responsive = TRUE,
      scrollZoom = TRUE
    )

  ma_interactive_file <- file.path(
    interactive_ma_dir,
    paste0(
      "MAplot_",
      comparison_name,
      "_interactive.html"
    )
  )

  save_interactive_widget(
    widget = p_ma_interactive,
    file = ma_interactive_file
  )

  # ==========================================================
  # 10.3 Échantillons de chaque groupe
  # ==========================================================

  numerator_samples <- colnames(
    normalized_counts
  )[
    sample_groups == numerator_group
  ]

  denominator_samples <- colnames(
    normalized_counts
  )[
    sample_groups == denominator_group
  ]

  if (
    length(numerator_samples) == 0 ||
      length(denominator_samples) == 0
  ) {

    stop(
      "Échantillons absents pour la comparaison ",
      comparison_name
    )
  }

  message(
    numerator_group,
    " : ",
    length(numerator_samples),
    " échantillons"
  )

  message(
    denominator_group,
    " : ",
    length(denominator_samples),
    " échantillons"
  )

  # ==========================================================
  # 10.4 Expression moyenne par groupe
  # ==========================================================

  numerator_mean <- rowMeans(
    normalized_counts[
      ,
      numerator_samples,
      drop = FALSE
    ]
  )

  denominator_mean <- rowMeans(
    normalized_counts[
      ,
      denominator_samples,
      drop = FALSE
    ]
  )

  scatter_df <- data.frame(
    gene = rownames(
      normalized_counts
    ),
    denominator_mean = denominator_mean,
    numerator_mean = numerator_mean,
    stringsAsFactors = FALSE
  )

  scatter_df$log_denominator <- log2(
    scatter_df$denominator_mean + 1
  )

  scatter_df$log_numerator <- log2(
    scatter_df$numerator_mean + 1
  )

  # ----------------------------------------------------------
  # Ajout des résultats différentiels
  # ----------------------------------------------------------

  result_match <- match(
    scatter_df$gene,
    res_df$gene
  )

  scatter_df$baseMean <- res_df$baseMean[
    result_match
  ]

  scatter_df$padj <- res_df$padj[
    result_match
  ]

  scatter_df$pvalue <- res_df$pvalue[
    result_match
  ]

  scatter_df$log2FoldChange <- res_df$log2FoldChange[
    result_match
  ]

  scatter_df$status <- "Non significatif"

  scatter_df$status[
    !is.na(scatter_df$padj) &
      scatter_df$padj < padj_threshold &
      !is.na(scatter_df$log2FoldChange) &
      scatter_df$log2FoldChange >= lfc_threshold
  ] <- paste0(
    "Plus élevé dans ",
    numerator_group
  )

  scatter_df$status[
    !is.na(scatter_df$padj) &
      scatter_df$padj < padj_threshold &
      !is.na(scatter_df$log2FoldChange) &
      scatter_df$log2FoldChange <= -lfc_threshold
  ] <- paste0(
    "Plus élevé dans ",
    denominator_group
  )

  scatter_df$status <- factor(
    scatter_df$status,
    levels = status_levels
  )

  # ----------------------------------------------------------
  # Texte interactif affiché au survol
  # ----------------------------------------------------------

  scatter_df$hover_text <- paste0(
    "<b>Gène : ",
    scatter_df$gene,
    "</b>",
    "<br>",
    denominator_group,
    " : ",
    format(
      round(
        scatter_df$denominator_mean,
        3
      ),
      big.mark = " ",
      scientific = FALSE
    ),
    "<br>",
    numerator_group,
    " : ",
    format(
      round(
        scatter_df$numerator_mean,
        3
      ),
      big.mark = " ",
      scientific = FALSE
    ),
    "<br>log2FC : ",
    round(
      scatter_df$log2FoldChange,
      3
    ),
    "<br>p-value : ",
    format(
      scatter_df$pvalue,
      scientific = TRUE,
      digits = 3
    ),
    "<br>padj : ",
    format(
      scatter_df$padj,
      scientific = TRUE,
      digits = 3
    ),
    "<br>Statut : ",
    scatter_df$status
  )

  # ==========================================================
  # 10.5 Corrélations
  # ==========================================================

  pearson_r <- cor(
    scatter_df$log_denominator,
    scatter_df$log_numerator,
    method = "pearson",
    use = "complete.obs"
  )

  spearman_rho <- cor(
    scatter_df$log_denominator,
    scatter_df$log_numerator,
    method = "spearman",
    use = "complete.obs"
  )

  # ==========================================================
  # 10.6 Scatterplot statique
  # ==========================================================

  scatter_file <- file.path(
    scatter_dir,
    paste0(
      "Scatterplot_",
      comparison_name,
      ".png"
    )
  )

  p_scatter <- ggplot(
    scatter_df,
    aes(
      x = log_denominator,
      y = log_numerator,
      color = status
    )
  ) +
    geom_point(
      alpha = 0.45,
      size = 1.2
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    coord_equal() +
    annotate(
      geom = "text",
      x = Inf,
      y = -Inf,
      hjust = 1.1,
      vjust = -0.7,
      label = paste0(
        "Pearson r = ",
        round(
          pearson_r,
          3
        ),
        "\nSpearman ρ = ",
        round(
          spearman_rho,
          3
        )
      ),
      size = 4.5
    ) +
    labs(
      title = gsub(
        "_",
        " ",
        comparison_name
      ),
      subtitle = paste0(
        "Moyenne des comptes normalisés par groupe",
        " — padj < ",
        padj_threshold,
        " et |log2FC| ≥ ",
        lfc_threshold
      ),
      x = paste0(
        denominator_group,
        "\nlog2(moyenne des comptes normalisés + 1)"
      ),
      y = paste0(
        numerator_group,
        "\nlog2(moyenne des comptes normalisés + 1)"
      ),
      color = "Résultat DESeq2"
    ) +
    theme_minimal(
      base_size = 13
    ) +
    theme(
      legend.position = "top",
      plot.title = element_text(
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      )
    )

  ggsave(
    filename = scatter_file,
    plot = p_scatter,
    width = 9,
    height = 8,
    dpi = 300
  )

  message(
    "Scatterplot statique sauvegardé : ",
    scatter_file
  )

  # ==========================================================
  # 10.7 Scatterplot interactif
  # ==========================================================

  interactive_axis_min <- min(
    c(
      scatter_df$log_denominator,
      scatter_df$log_numerator
    ),
    na.rm = TRUE
  )

  interactive_axis_max <- max(
    c(
      scatter_df$log_denominator,
      scatter_df$log_numerator
    ),
    na.rm = TRUE
  )

  p_scatter_interactive <- plot_ly(
    data = scatter_df,
    x = ~log_denominator,
    y = ~log_numerator,
    color = ~status,
    colors = "Set1",
    type = "scatter",
    mode = "markers",
    text = ~hover_text,
    hoverinfo = "text",
    marker = list(
      size = 5,
      opacity = 0.55
    )
  )

  p_scatter_interactive <- p_scatter_interactive |>
    layout(
      title = list(
        text = paste0(
          gsub(
            "_",
            " ",
            comparison_name
          ),
          "<br>",
          "<sup>",
          "Pearson r = ",
          round(
            pearson_r,
            3
          ),
          " — Spearman ρ = ",
          round(
            spearman_rho,
            3
          ),
          "</sup>"
        ),
        x = 0.5
      ),
      xaxis = list(
        title = paste0(
          denominator_group,
          "<br>",
          "log2(moyenne des comptes normalisés + 1)"
        ),
        range = c(
          interactive_axis_min,
          interactive_axis_max
        ),
        scaleanchor = "y",
        scaleratio = 1,
        zeroline = FALSE
      ),
      yaxis = list(
        title = paste0(
          numerator_group,
          "<br>",
          "log2(moyenne des comptes normalisés + 1)"
        ),
        range = c(
          interactive_axis_min,
          interactive_axis_max
        ),
        zeroline = FALSE
      ),
      legend = list(
        title = list(
          text = "Résultat DESeq2"
        ),
        orientation = "h",
        x = 0,
        y = 1.15
      ),
      shapes = list(
        list(
          type = "line",
          x0 = interactive_axis_min,
          x1 = interactive_axis_max,
          y0 = interactive_axis_min,
          y1 = interactive_axis_max,
          line = list(
            dash = "dash",
            width = 1
          )
        )
      ),
      hovermode = "closest"
    ) |>
    config(
      displaylogo = FALSE,
      responsive = TRUE,
      scrollZoom = TRUE
    )

  scatter_interactive_file <- file.path(
    interactive_scatter_dir,
    paste0(
      "Scatterplot_",
      comparison_name,
      "_interactive.html"
    )
  )

  save_interactive_widget(
    widget = p_scatter_interactive,
    file = scatter_interactive_file
  )

  # ==========================================================
  # 10.8 Sauvegarde des données du scatterplot
  # ==========================================================

  scatter_csv <- file.path(
    scatter_dir,
    paste0(
      "Scatterplot_data_",
      comparison_name,
      ".csv"
    )
  )

  write.csv(
    scatter_df,
    scatter_csv,
    row.names = FALSE
  )

  message(
    "Données du scatterplot sauvegardées : ",
    scatter_csv
  )

  # ==========================================================
  # 10.9 Résumé
  # ==========================================================

  significant_up <- sum(
    !is.na(res_df$padj) &
      res_df$padj < padj_threshold &
      !is.na(res_df$log2FoldChange) &
      res_df$log2FoldChange >= lfc_threshold
  )

  significant_down <- sum(
    !is.na(res_df$padj) &
      res_df$padj < padj_threshold &
      !is.na(res_df$log2FoldChange) &
      res_df$log2FoldChange <= -lfc_threshold
  )

  summary_list[[comparison_name]] <- data.frame(
    comparison = comparison_name,
    numerator = numerator_group,
    denominator = denominator_group,
    numerator_samples = length(
      numerator_samples
    ),
    denominator_samples = length(
      denominator_samples
    ),
    pearson_r = pearson_r,
    spearman_rho = spearman_rho,
    significant_higher_numerator = significant_up,
    significant_higher_denominator = significant_down,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 11. Tableau résumé
# ============================================================

comparison_summary <- do.call(
  rbind,
  summary_list
)

rownames(
  comparison_summary
) <- NULL

summary_file <- file.path(
  script_dir,
  "Comparison_summary.csv"
)

write.csv(
  comparison_summary,
  summary_file,
  row.names = FALSE
)

message(
  "Résumé sauvegardé : ",
  summary_file
)

print(
  comparison_summary
)

message(
  "Toutes les analyses sont terminées."
)

message(
  "Figures interactives disponibles dans : ",
  interactive_dir
)
