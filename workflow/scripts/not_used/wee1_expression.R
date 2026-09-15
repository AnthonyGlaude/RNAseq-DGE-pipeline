#!/usr/bin/env Rscript

# ============================================================
# 0. Packages
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(ggrepel)
  library(patchwork)
})


# ============================================================
# 1. Chemins
# ============================================================

base_dir <- paste0(
  "/home/glaudea/scratch/glaudea/test_souris_rnaseq/",
  "RNA_seq-analysis/workflow/results_avec_wee1as"
)

dds_file <- file.path(
  base_dir,
  "deseq2",
  "dds.rds"
)

plot_dir <- file.path(
  base_dir,
  "wee1_expression"
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Paramètres globaux
# ============================================================

group_levels <- c(
  "Control_WT",
  "Control_KO",
  "HFpEF_WT",
  "HFpEF_KO"
)

group_colors <- c(
  "Control_WT" = "#0072B2",
  "Control_KO" = "#56B4E9",
  "HFpEF_WT"   = "#D55E00",
  "HFpEF_KO"   = "#E69F00"
)


# ============================================================
# 3. Charger le DESeqDataSet
# ============================================================

if (!file.exists(dds_file)) {
  stop(
    "Fichier DDS introuvable : ",
    dds_file
  )
}

dds <- readRDS(dds_file)

message("")
message("========================================")
message("DDS chargé")
message("========================================")

message(
  "Nombre de gènes : ",
  nrow(dds)
)

message(
  "Nombre d'échantillons : ",
  ncol(dds)
)

metadata <- as.data.frame(
  colData(dds)
)

message("")
message("Métadonnées du DDS :")

print(metadata)


# ============================================================
# 4. Construire les 4 groupes
# ============================================================

group_column <- NULL

for (column_name in colnames(metadata)) {

  column_values <- unique(
    as.character(
      metadata[[column_name]]
    )
  )

  if (all(group_levels %in% column_values)) {

    group_column <- column_name

    break
  }
}


if (!is.null(group_column)) {

  message(
    "Les 4 groupes ont été trouvés dans la colonne : ",
    group_column
  )

  colData(dds)$group <- factor(
    as.character(
      colData(dds)[[group_column]]
    ),
    levels = group_levels
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
    "Création de 'group' à partir de condition + genotype."
  )


  # ----------------------------------------------------------
  # Condition
  # ----------------------------------------------------------

  condition_clean <- trimws(
    as.character(
      colData(dds)$condition
    )
  )

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
      "hfpef|patho",
      condition_clean,
      ignore.case = TRUE
    )
  ] <- "HFpEF"


  # ----------------------------------------------------------
  # Génotype
  # ----------------------------------------------------------

  genotype_clean <- trimws(
    as.character(
      colData(dds)$genotype
    )
  )

  genotype_clean[
    grepl(
      "^wt$|wild",
      genotype_clean,
      ignore.case = TRUE
    )
  ] <- "WT"

  genotype_clean[
    grepl(
      "^ko$|knockout|knock-out",
      genotype_clean,
      ignore.case = TRUE
    )
  ] <- "KO"


  # ----------------------------------------------------------
  # Groupe
  # ----------------------------------------------------------

  colData(dds)$group <- factor(
    paste(
      condition_clean,
      genotype_clean,
      sep = "_"
    ),
    levels = group_levels
  )

} else {

  stop(
    paste0(
      "Impossible de reconstruire les 4 groupes.\n",
      "Il faut une colonne contenant directement les groupes ",
      "ou les colonnes condition + genotype."
    )
  )
}


message("")
message("Répartition des groupes :")

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

  stop(
    "Au moins un échantillon n'a pas été assigné à un groupe."
  )
}


# ============================================================
# 5. Construire condition et genotype propres
# ============================================================

colData(dds)$condition_analysis <- factor(

  ifelse(
    colData(dds)$group %in%
      c(
        "Control_WT",
        "Control_KO"
      ),
    "Control",
    "HFpEF"
  ),

  levels = c(
    "Control",
    "HFpEF"
  )
)


colData(dds)$genotype_analysis <- factor(

  ifelse(
    colData(dds)$group %in%
      c(
        "Control_WT",
        "HFpEF_WT"
      ),
    "WT",
    "KO"
  ),

  levels = c(
    "WT",
    "KO"
  )
)


# ============================================================
# 6. DESeq2 avec interaction condition × genotype
# ============================================================

design(dds) <-
  ~ condition_analysis * genotype_analysis

dds <- DESeq(
  dds,
  quiet = TRUE
)


message("")
message("========================================")
message("Coefficients DESeq2")
message("========================================")

print(
  resultsNames(dds)
)


# ============================================================
# 7. Identifier les coefficients
# ============================================================

coef_names <- resultsNames(dds)


condition_coef <- grep(
  "^condition_analysis_.*_vs_.*$",
  coef_names,
  value = TRUE
)

genotype_coef <- grep(
  "^genotype_analysis_.*_vs_.*$",
  coef_names,
  value = TRUE
)

interaction_coef <- grep(
  "condition_analysis.*genotype_analysis",
  coef_names,
  value = TRUE
)


if (length(condition_coef) != 1) {

  stop(
    "Coefficient condition introuvable ou ambigu : ",
    paste(
      condition_coef,
      collapse = ", "
    )
  )
}


if (length(genotype_coef) != 1) {

  stop(
    "Coefficient genotype introuvable ou ambigu : ",
    paste(
      genotype_coef,
      collapse = ", "
    )
  )
}


if (length(interaction_coef) != 1) {

  stop(
    "Coefficient interaction introuvable ou ambigu : ",
    paste(
      interaction_coef,
      collapse = ", "
    )
  )
}


message("")
message(
  "Condition   : ",
  condition_coef
)

message(
  "Genotype    : ",
  genotype_coef
)

message(
  "Interaction : ",
  interaction_coef
)


# ============================================================
# 8. Résultats DESeq2 des 4 comparaisons
# ============================================================


# ------------------------------------------------------------
# 8.1 Control KO vs Control WT
# ------------------------------------------------------------

res_control_ko_vs_wt <- results(
  dds,
  name = genotype_coef
)


# ------------------------------------------------------------
# 8.2 HFpEF WT vs Control WT
# ------------------------------------------------------------

res_hfpef_wt_vs_control_wt <- results(
  dds,
  name = condition_coef
)


# ------------------------------------------------------------
# 8.3 HFpEF KO vs HFpEF WT
# genotype + interaction
# ------------------------------------------------------------

res_hfpef_ko_vs_wt <- results(
  dds,
  contrast = list(
    c(
      genotype_coef,
      interaction_coef
    )
  )
)


# ------------------------------------------------------------
# 8.4 HFpEF KO vs Control KO
# condition + interaction
# ------------------------------------------------------------

res_hfpef_ko_vs_control_ko <- results(
  dds,
  contrast = list(
    c(
      condition_coef,
      interaction_coef
    )
  )
)


# ============================================================
# 9. Normalized counts
# ============================================================

norm_counts <- counts(
  dds,
  normalized = TRUE
)


# ============================================================
# 10. Identifier Wee1-AS
# ============================================================

gene_ids <- rownames(dds)


wee1as_candidates <- gene_ids[
  grepl(
    "wee1.?as",
    gene_ids,
    ignore.case = TRUE
  )
]


message("")
message("========================================")
message("Candidats Wee1-AS")
message("========================================")

print(
  wee1as_candidates
)


if (length(wee1as_candidates) == 0) {

  stop(
    paste0(
      "Wee1-AS n'a pas été trouvé dans le DDS.\n",
      "Vérifie le gene_id du GTF custom."
    )
  )
}


if (length(wee1as_candidates) > 1) {

  warning(
    "Plusieurs candidats Wee1-AS trouvés. ",
    "Le premier sera utilisé."
  )
}


wee1as_id <- wee1as_candidates[1]


message(
  "Wee1-AS utilisé : ",
  wee1as_id
)


# ============================================================
# 11. Identifier Wee1 canonique
# ============================================================

dds_ids_no_version <- sub(
  "\\..*$",
  "",
  gene_ids
)


ensembl_ids <- unique(

  dds_ids_no_version[
    grepl(
      "^ENSMUSG",
      dds_ids_no_version
    )
  ]
)


annotation <- AnnotationDbi::select(

  org.Mm.eg.db,

  keys = ensembl_ids,

  keytype = "ENSEMBL",

  columns = "SYMBOL"

) %>%

  dplyr::filter(
    !is.na(SYMBOL)
  ) %>%

  dplyr::distinct(
    ENSEMBL,
    .keep_all = TRUE
  )


wee1_annotation <- annotation %>%

  dplyr::filter(
    SYMBOL == "Wee1"
  )


if (nrow(wee1_annotation) == 0) {

  stop(
    "Impossible de trouver Wee1 dans org.Mm.eg.db."
  )
}


wee1_ensembl <-
  wee1_annotation$ENSEMBL[1]


wee1_match <- which(
  dds_ids_no_version ==
    wee1_ensembl
)


if (length(wee1_match) == 0) {

  stop(
    paste0(
      "Wee1 (",
      wee1_ensembl,
      ") est absent du DDS."
    )
  )
}


wee1_id <- gene_ids[
  wee1_match[1]
]


message("")
message("========================================")
message("Identifiants retenus")
message("========================================")

message(
  "Wee1-AS : ",
  wee1as_id
)

message(
  "Wee1    : ",
  wee1_id
)


# ============================================================
# 12. Fonction récupération statistiques DESeq2
# ============================================================

get_gene_stats <- function(
  gene_id
) {

  data.frame(

    comparison = c(
      "Control WT vs KO",
      "HFpEF WT vs KO",
      "WT Control vs HFpEF",
      "KO Control vs HFpEF"
    ),


    padj = c(

      res_control_ko_vs_wt[
        gene_id,
        "padj"
      ],

      res_hfpef_ko_vs_wt[
        gene_id,
        "padj"
      ],

      res_hfpef_wt_vs_control_wt[
        gene_id,
        "padj"
      ],

      res_hfpef_ko_vs_control_ko[
        gene_id,
        "padj"
      ]
    ),


    log2FC = c(

      res_control_ko_vs_wt[
        gene_id,
        "log2FoldChange"
      ],

      res_hfpef_ko_vs_wt[
        gene_id,
        "log2FoldChange"
      ],

      res_hfpef_wt_vs_control_wt[
        gene_id,
        "log2FoldChange"
      ],

      res_hfpef_ko_vs_control_ko[
        gene_id,
        "log2FoldChange"
      ]
    )
  )
}


# ============================================================
# 13. Fonction étoiles
# ============================================================

p_to_signif <- function(p) {

  if (is.na(p)) {
    return("NA")
  }

  if (p < 0.0001) {
    return("****")
  }

  if (p < 0.001) {
    return("***")
  }

  if (p < 0.01) {
    return("**")
  }

  if (p < 0.05) {
    return("*")
  }

  return("ns")
}


# ============================================================
# 14. Tableau d'expression
# ============================================================

expression_df <- data.frame(

  sample = colnames(dds),

  group = factor(
    as.character(
      colData(dds)$group
    ),
    levels = group_levels
  ),

  Wee1_AS = as.numeric(
    norm_counts[
      wee1as_id,
    ]
  ),

  Wee1 = as.numeric(
    norm_counts[
      wee1_id,
    ]
  )
)


# Numéro de souris
expression_df$sample_label <- sub(
  "^.*_",
  "",
  expression_df$sample
)


message("")
message("========================================")
message("Expression individuelle")
message("========================================")

print(
  expression_df
)


# ============================================================
# 15. Format long
# ============================================================

expression_long <- expression_df %>%

  tidyr::pivot_longer(

    cols = c(
      Wee1_AS,
      Wee1
    ),

    names_to = "gene",

    values_to = "normalized_count"

  ) %>%

  dplyr::mutate(

    gene = dplyr::recode(
      gene,
      "Wee1_AS" = "Wee1-AS",
      "Wee1" = "Wee1"
    )
  )


# ============================================================
# 16. Moyennes par groupe
# ============================================================

summary_df <- expression_long %>%

  dplyr::group_by(
    gene,
    group
  ) %>%

  dplyr::summarise(

    mean_normalized_count = mean(
      normalized_count,
      na.rm = TRUE
    ),

    sd_normalized_count = sd(
      normalized_count,
      na.rm = TRUE
    ),

    n = dplyr::n(),

    .groups = "drop"
  )


message("")
message("========================================")
message("Moyennes par groupe")
message("========================================")

print(
  summary_df
)


# ============================================================
# 17. Statistiques Wee1-AS
# ============================================================

stats_wee1as <- get_gene_stats(
  wee1as_id
)


stats_wee1as$significance <- vapply(

  stats_wee1as$padj,

  p_to_signif,

  character(1)
)


message("")
message("========================================")
message("DESeq2 : Wee1-AS")
message("========================================")

print(
  stats_wee1as
)


# ============================================================
# 18. Statistiques Wee1
# ============================================================

stats_wee1 <- get_gene_stats(
  wee1_id
)


stats_wee1$significance <- vapply(

  stats_wee1$padj,

  p_to_signif,

  character(1)
)


message("")
message("========================================")
message("DESeq2 : Wee1")
message("========================================")

print(
  stats_wee1
)


# ============================================================
# 19. ÉCHELLE Y COMMUNE POUR WEE1-AS ET WEE1
# ============================================================
#
# Très important :
# les deux graphiques utilisent exactement les mêmes limites Y.
#
# ============================================================

shared_y_max_data <- max(
  expression_long$normalized_count,
  na.rm = TRUE
)

shared_y_min <- 0

shared_y_range <- shared_y_max_data -
  shared_y_min


if (
  !is.finite(shared_y_range) ||
  shared_y_range == 0
) {

  shared_y_range <- 1
}


# ------------------------------------------------------------
# Positions communes des crochets
# ------------------------------------------------------------

shared_y_1 <-
  shared_y_max_data +
  0.12 * shared_y_range

shared_y_2 <-
  shared_y_max_data +
  0.27 * shared_y_range

shared_y_3 <-
  shared_y_max_data +
  0.42 * shared_y_range


shared_bracket_height <-
  0.025 * shared_y_range


# Limite supérieure finale de la figure
shared_plot_y_max <-
  shared_y_max_data +
  0.55 * shared_y_range


message("")
message("========================================")
message("Échelle Y commune")
message("========================================")

message(
  "Maximum observé : ",
  round(
    shared_y_max_data,
    2
  )
)

message(
  "Limite Y commune : 0 - ",
  round(
    shared_plot_y_max,
    2
  )
)


# ============================================================
# 20. Fonction de boxplot
# ============================================================

make_boxplot <- function(
  gene_name,
  stats_gene
) {


  # ----------------------------------------------------------
  # Données du gène
  # ----------------------------------------------------------

  gene_df <- expression_long %>%

    dplyr::filter(
      gene == gene_name
    )


  # ----------------------------------------------------------
  # Significativité
  # ----------------------------------------------------------

  sig_control_wt_ko <-

    stats_gene$significance[
      stats_gene$comparison ==
        "Control WT vs KO"
    ]


  sig_hfpef_wt_ko <-

    stats_gene$significance[
      stats_gene$comparison ==
        "HFpEF WT vs KO"
    ]


  sig_wt_control_hf <-

    stats_gene$significance[
      stats_gene$comparison ==
        "WT Control vs HFpEF"
    ]


  sig_ko_control_hf <-

    stats_gene$significance[
      stats_gene$comparison ==
        "KO Control vs HFpEF"
    ]


  # ----------------------------------------------------------
  # Même position des crochets pour les deux gènes
  # ----------------------------------------------------------

  y_1 <- shared_y_1
  y_2 <- shared_y_2
  y_3 <- shared_y_3

  bracket_height <-
    shared_bracket_height


  # ----------------------------------------------------------
  # Jitter reproductible
  # ----------------------------------------------------------

#  jitter_position <- position_jitter(
#    width = 0.07,
#    height = 0,
#    seed = 123
#  )


  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------

  p_box <- ggplot(

    gene_df,

    aes(
      x = group,
      y = normalized_count,
      fill = group,
      color = group
    )

  ) +


    # ========================================================
    # Boxplot
    # ========================================================

    geom_boxplot(
      width = 0.45,
      alpha = 0.20,
      outlier.shape = NA,
      linewidth = 0.7
    ) +


    # ========================================================
    # Points individuels
    # ========================================================

    geom_point(
    #  position = jitter_position,
      size = 2.7
    ) +


    # ========================================================
    # Numéros des souris
    # ========================================================

    geom_text_repel(

      aes(
        label = sample_label
      ),

     # position = jitter_position,

      size = 3.2,

      direction = "y",

      box.padding = 0.25,

      point.padding = 0.20,
      segment.color = NA,
      min.segment.length = 0,

      max.overlaps = Inf,

      show.legend = FALSE
    ) +


    # ========================================================
    # 1. Control WT vs Control KO
    # ========================================================

    annotate(
      "segment",
      x = 1,
      xend = 2,
      y = y_1,
      yend = y_1,
      color = "black",
      linewidth = 0.6
    ) +

    annotate(
      "segment",
      x = 1,
      xend = 1,
      y = y_1,
      yend = y_1 - bracket_height,
      color = "black",
      linewidth = 0.6
    ) +

    annotate(
      "segment",
      x = 2,
      xend = 2,
      y = y_1,
      yend = y_1 - bracket_height,
      color = "black",
      linewidth = 0.6
    ) +

    annotate(
      "text",
      x = 1.5,
      y = y_1 +
        0.025 * shared_y_range,
      label = sig_control_wt_ko,
      color = "black",
      size = 4.2
    ) +


    # ========================================================
    # 2. HFpEF WT vs HFpEF KO
    # ========================================================

    annotate(
      "segment",
      x = 3,
      xend = 4,
      y = y_1,
      yend = y_1,
      color = "black",
      linewidth = 0.6
    ) +

    annotate(
      "segment",
      x = 3,
      xend = 3,
      y = y_1,
      yend = y_1 - bracket_height,
      color = "black",
      linewidth = 0.6
    ) +

    annotate(
      "segment",
      x = 4,
      xend = 4,
      y = y_1,
      yend = y_1 - bracket_height,
      color = "black",
      linewidth = 0.6
    ) +

    annotate(
      "text",
      x = 3.5,
      y = y_1 +
        0.025 * shared_y_range,
      label = sig_hfpef_wt_ko,
      color = "black",
      size = 4.2
    ) +


    # ========================================================
    # 3. Control WT vs HFpEF WT
    # ========================================================

    annotate(
      "segment",
      x = 1,
      xend = 3,
      y = y_2,
      yend = y_2,
      color = "black",
      linewidth = 0.6
    ) +

    annotate(
      "segment",
      x = 1,
      xend = 1,
      y = y_2,
      yend = y_2 - bracket_height,
      color = "black",
      linewidth = 0.6
    ) +

    annotate(
      "segment",
      x = 3,
      xend = 3,
      y = y_2,
      yend = y_2 - bracket_height,
      color = "black",
      linewidth = 0.6
    ) +

    annotate(
      "text",
      x = 2,
      y = y_2 +
        0.025 * shared_y_range,
      label = sig_wt_control_hf,
      color = "black",
      size = 4.2
    ) +


    # ========================================================
    # 4. Control KO vs HFpEF KO
    # ========================================================

    annotate(
      "segment",
      x = 2,
      xend = 4,
      y = y_3,
      yend = y_3,
      color = "black",
      linewidth = 0.6
    ) +

    annotate(
      "segment",
      x = 2,
      xend = 2,
      y = y_3,
      yend = y_3 - bracket_height,
      color = "black",
      linewidth = 0.6
    ) +

    annotate(
      "segment",
      x = 4,
      xend = 4,
      y = y_3,
      yend = y_3 - bracket_height,
      color = "black",
      linewidth = 0.6
    ) +

    annotate(
      "text",
      x = 3,
      y = y_3 +
        0.025 * shared_y_range,
      label = sig_ko_control_hf,
      color = "black",
      size = 4.2
    ) +


    # ========================================================
    # Couleurs
    # ========================================================

    scale_color_manual(
      values = group_colors,
      breaks = group_levels,
      drop = FALSE
    ) +

    scale_fill_manual(
      values = group_colors,
      breaks = group_levels,
      drop = FALSE
    ) +


    # ========================================================
    # Labels axe X
    # ========================================================

    scale_x_discrete(

      labels = c(
        "Control_WT" = "Control WT",
        "Control_KO" = "Control KO",
        "HFpEF_WT"   = "HFpEF WT",
        "HFpEF_KO"   = "HFpEF KO"
      )
    ) +


    # ========================================================
    # ÉCHELLE Y COMMUNE
    # ========================================================

    scale_y_continuous(

      limits = c(
        shared_y_min,
        shared_plot_y_max
      ),

      expand = expansion(
        mult = c(
          0.02,
          0.02
        )
      )
    ) +


    # ========================================================
    # Labels
    # ========================================================

    labs(

      title = paste0(
        "Expression of ",
        gene_name
      ),

      x = NULL,

      y = "DESeq2 normalized counts",

      caption =
        "DESeq2 Wald test; BH-adjusted p-values"
    ) +


    # ========================================================
    # Thème
    # ========================================================

    theme_classic(
      base_size = 13
    ) +

    theme(

      legend.position = "none",

      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      ),

      plot.caption = element_text(
        hjust = 0,
        size = 9,
        margin = margin(
          t = 10
        )
      ),

      axis.title.y = element_text(
        size = 12
      ),

      axis.text.x = element_text(
        angle = 20,
        hjust = 1,
        size = 11
      ),

      axis.text.y = element_text(
        size = 10
      )
    )


  return(
    p_box
  )
}


# ============================================================
# 21. Générer Wee1-AS
# ============================================================

p_wee1as <- make_boxplot(

  gene_name = "Wee1-AS",

  stats_gene = stats_wee1as
)


print(
  p_wee1as
)


wee1as_file <- file.path(
  plot_dir,
  "Wee1-AS_expression_boxplot_4groups.png"
)


ggsave(

  filename = wee1as_file,

  plot = p_wee1as,

  width = 8,

  height = 6.5,

  dpi = 300,

  bg = "white"
)


# ============================================================
# 22. Générer Wee1
# ============================================================

p_wee1 <- make_boxplot(

  gene_name = "Wee1",

  stats_gene = stats_wee1
)


print(
  p_wee1
)


wee1_file <- file.path(
  plot_dir,
  "Wee1_expression_boxplot_4groups.png"
)


ggsave(

  filename = wee1_file,

  plot = p_wee1,

  width = 8,

  height = 6.5,

  dpi = 300,

  bg = "white"
)


# ============================================================
# 23. Figure combinée
# ============================================================

p_combined <- (

  p_wee1as +

    p_wee1 +

    patchwork::plot_layout(
      ncol = 2
    )

) +

  patchwork::plot_annotation(

    title =
      " "
  )


combined_file <- file.path(

  plot_dir,

  "Wee1-AS_Wee1_expression_boxplots_combined_sameY.png"
)


ggsave(

  filename = combined_file,

  plot = p_combined,

  width = 16,

  height = 6.5,

  dpi = 300,

  bg = "white"
)

# ============================================================
# 24. Corrélation Wee1-AS vs Wee1
# ============================================================

message("")
message("========================================")
message("CORRÉLATION WEE1-AS vs WEE1")
message("========================================")

# ------------------------------------------------------------
# Ajouter condition et génotype au tableau
# ------------------------------------------------------------

expression_df <- expression_df %>%
  mutate(
    condition = ifelse(
      group %in% c("Control_WT", "Control_KO"),
      "Control",
      "HFpEF"
    ),
    genotype = ifelse(
      group %in% c("Control_WT", "HFpEF_WT"),
      "WT",
      "KO"
    ),
    condition = factor(
      condition,
      levels = c("Control", "HFpEF")
    ),
    genotype = factor(
      genotype,
      levels = c("WT", "KO")
    )
  )


# ============================================================
# 24.1 Corrélation globale
# ============================================================

cor_global_pearson <- cor.test(
  expression_df$Wee1_AS,
  expression_df$Wee1,
  method = "pearson"
)

cor_global_spearman <- cor.test(
  expression_df$Wee1_AS,
  expression_df$Wee1,
  method = "spearman",
  exact = FALSE
)

message("")
message("GLOBAL")
message(
  "Pearson r = ",
  round(unname(cor_global_pearson$estimate), 3),
  " ; p = ",
  signif(cor_global_pearson$p.value, 3)
)

message(
  "Spearman rho = ",
  round(unname(cor_global_spearman$estimate), 3),
  " ; p = ",
  signif(cor_global_spearman$p.value, 3)
)


# ============================================================
# 24.2 Corrélation séparée Control / HFpEF
# ============================================================

cor_by_condition <- expression_df %>%
  group_by(condition) %>%
  group_modify(
    ~ {
      pearson <- cor.test(
        .x$Wee1_AS,
        .x$Wee1,
        method = "pearson"
      )

      spearman <- cor.test(
        .x$Wee1_AS,
        .x$Wee1,
        method = "spearman",
        exact = FALSE
      )

      tibble(
        n = nrow(.x),
        pearson_r = unname(pearson$estimate),
        pearson_p = pearson$p.value,
        spearman_rho = unname(spearman$estimate),
        spearman_p = spearman$p.value
      )
    }
  )

message("")
message("PAR CONDITION")
print(cor_by_condition)


# ============================================================
# 24.3 Corrélation ajustée pour les 4 groupes
# ============================================================
#
# On retire de chaque variable la différence moyenne attribuable
# au groupe expérimental.
#
# On corrèle ensuite les résidus :
# = variation Wee1-AS et Wee1 ENTRE individus au-delà du groupe.
#
# ============================================================

model_wee1as <- lm(
  Wee1_AS ~ group,
  data = expression_df
)

model_wee1 <- lm(
  Wee1 ~ group,
  data = expression_df
)

expression_df$Wee1_AS_residual <- residuals(model_wee1as)
expression_df$Wee1_residual <- residuals(model_wee1)

cor_adjusted_pearson <- cor.test(
  expression_df$Wee1_AS_residual,
  expression_df$Wee1_residual,
  method = "pearson"
)

cor_adjusted_spearman <- cor.test(
  expression_df$Wee1_AS_residual,
  expression_df$Wee1_residual,
  method = "spearman",
  exact = FALSE
)

message("")
message("AJUSTÉ POUR LES 4 GROUPES")
message(
  "Pearson r = ",
  round(unname(cor_adjusted_pearson$estimate), 3),
  " ; p = ",
  signif(cor_adjusted_pearson$p.value, 3)
)

message(
  "Spearman rho = ",
  round(unname(cor_adjusted_spearman$estimate), 3),
  " ; p = ",
  signif(cor_adjusted_spearman$p.value, 3)
)


# ============================================================
# 24.4 Figure : corrélation brute
# ============================================================

global_label <- paste0(
  "Pearson r = ",
  round(unname(cor_global_pearson$estimate), 2),
  "\nP = ",
  format.pval(
    cor_global_pearson$p.value,
    digits = 2,
    eps = 0.0001
  )
)

p_correlation_global <- ggplot(
  expression_df,
  aes(
    x = Wee1_AS,
    y = Wee1,
    color = group
  )
) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    fill = "grey80",
    linewidth = 0.8
  ) +
  geom_point(
    size = 3.5
  ) +
  geom_text_repel(
    aes(label = sample_label),
    size = 3.2,
    box.padding = 0.30,
    point.padding = 0.20,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = global_label,
    hjust = -0.10,
    vjust = 1.15,
    size = 4
  ) +
  scale_color_manual(
    values = group_colors,
    breaks = group_levels,
    drop = FALSE
  ) +
  labs(
    title = "Corrélation brute entre Wee1-AS et Wee1",
    x = "Wee1-AS — normalized counts",
    y = "Wee1 — normalized counts",
    color = "Groupe"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.position = "top"
  )


# ============================================================
# 24.5 Figure : corrélation ajustée
# ============================================================

adjusted_label <- paste0(
  "Pearson r = ",
  round(unname(cor_adjusted_pearson$estimate), 2),
  "\nP = ",
  format.pval(
    cor_adjusted_pearson$p.value,
    digits = 2,
    eps = 0.0001
  )
)

p_correlation_adjusted <- ggplot(
  expression_df,
  aes(
    x = Wee1_AS_residual,
    y = Wee1_residual,
    color = group
  )
) +
  geom_hline(
    yintercept = 0,
    color = "grey75",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    color = "grey75",
    linewidth = 0.4
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    fill = "grey80",
    linewidth = 0.8
  ) +
  geom_point(
    size = 3.5
  ) +
  geom_text_repel(
    aes(label = sample_label),
    size = 3.2,
    box.padding = 0.30,
    point.padding = 0.20,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = adjusted_label,
    hjust = -0.10,
    vjust = 1.15,
    size = 4
  ) +
  scale_color_manual(
    values = group_colors,
    breaks = group_levels,
    drop = FALSE
  ) +
  labs(
    title = "Corrélation après ajustement pour le groupe",
    x = "Wee1-AS — résidus",
    y = "Wee1 — résidus",
    color = "Groupe"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.position = "top"
  )


# ============================================================
# 24.6 Figure combinée
# ============================================================

p_correlation_combined <-
  p_correlation_global +
  p_correlation_adjusted +
  patchwork::plot_layout(ncol = 2)

print(p_correlation_combined)


# ============================================================
# 24.7 Sauvegarder
# ============================================================

ggsave(
  filename = file.path(
    plot_dir,
    "Wee1-AS_Wee1_correlation_global_vs_adjusted.png"
  ),
  plot = p_correlation_combined,
  width = 14,
  height = 6.5,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 24.8 Export statistiques
# ============================================================

correlation_stats <- bind_rows(

  data.frame(
    analysis = "Global",
    method = c("Pearson", "Spearman"),
    correlation = c(
      unname(cor_global_pearson$estimate),
      unname(cor_global_spearman$estimate)
    ),
    pvalue = c(
      cor_global_pearson$p.value,
      cor_global_spearman$p.value
    )
  ),

  data.frame(
    analysis = "Adjusted_for_4_groups",
    method = c("Pearson", "Spearman"),
    correlation = c(
      unname(cor_adjusted_pearson$estimate),
      unname(cor_adjusted_spearman$estimate)
    ),
    pvalue = c(
      cor_adjusted_pearson$p.value,
      cor_adjusted_spearman$p.value
    )
  )

)

write.csv(
  correlation_stats,
  file.path(
    plot_dir,
    "Wee1-AS_Wee1_correlation_statistics.csv"
  ),
  row.names = FALSE
)

write.csv(
  cor_by_condition,
  file.path(
    plot_dir,
    "Wee1-AS_Wee1_correlation_by_condition.csv"
  ),
  row.names = FALSE
)
