#!/usr/bin/env Rscript


# ============================================================
# 0. LIBRARIES
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(enrichplot)
  library(ggtext)
})


# ============================================================
# 1. PATHS
# ============================================================

stats_file <- paste0(
  "/home/glaudea/scratch/glaudea/test_souris_rnaseq/",
  "RNA_seq-analysis/workflow/resultat_article/results/",
  "deseq2_WT_only/",
  "HFpEF_WT-Control_WT_DESeq2_gene.csv"
)


# IMPORTANT:
# ONE common GO source for BP + CC + MF
gmt_file <- paste0(
  "/home/glaudea/scratch/glaudea/test_souris_rnaseq/",
  "RNA_seq-analysis/workflow/resultat_article/results/",
  "deseq2_WT_only/",
  "Mouse_GO_AllPathways_noPFOCRno_GO_iea_August_10_2026_symbol.gmt"
)


out_dir <- paste0(
  "/home/glaudea/scratch/glaudea/test_souris_rnaseq/",
  "RNA_seq-analysis/workflow/resultat_article/results/",
  "GSEA_WT_only"
)


dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. PARAMETERS
# ============================================================

gsea_fdr_cutoff <- 0.25

min_gs_size <- 15
max_gs_size <- 500

n_top_up <- 10
n_top_down <- 10
n_top_cc_negative <- 15
n_top_kegg_up <- 10
n_top_kegg_down <- 10


# ============================================================
# 3. CHECK FILES
# ============================================================

if (!file.exists(stats_file)) {
  stop(
    paste0(
      "DESeq2 file not found:\n",
      stats_file
    )
  )
}


if (!file.exists(gmt_file)) {
  stop(
    paste0(
      "GMT file not found:\n",
      gmt_file
    )
  )
}


# ============================================================
# 4. READ DESEQ2 RESULTS
# ============================================================

df <- read_csv(
  stats_file,
  show_col_types = FALSE
)


cat("\nColumns in DESeq2 file:\n")
print(names(df))


# ============================================================
# 5. DETECT ENSEMBL COLUMN
# ============================================================

candidate_gene_columns <- c(
  "gene",
  "Gene",
  "ENSEMBL",
  "ensembl",
  "ensembl_gene_id",
  "gene_id",
  "...1"
)


gene_col <- candidate_gene_columns[
  candidate_gene_columns %in% names(df)
][1]


if (is.na(gene_col)) {

  ensembl_match <- vapply(
    df,
    function(x) {

      x <- as.character(x)

      any(
        grepl(
          "^ENSMUSG",
          x
        ),
        na.rm = TRUE
      )
    },
    logical(1)
  )


  matching_columns <- names(df)[
    ensembl_match
  ]


  if (length(matching_columns) == 0) {

    stop(
      paste0(
        "Could not identify an Ensembl gene column.\n",
        "Available columns:\n",
        paste(
          names(df),
          collapse = ", "
        )
      )
    )
  }


  gene_col <- matching_columns[1]
}


cat(
  "\nEnsembl column detected: ",
  gene_col,
  "\n",
  sep = ""
)


# ============================================================
# 6. CHECK REQUIRED DESEQ2 COLUMNS
# ============================================================

required_columns <- c(
  "log2FoldChange",
  "stat"
)


missing_columns <- setdiff(
  required_columns,
  names(df)
)


if (length(missing_columns) > 0) {

  stop(
    paste0(
      "Missing DESeq2 columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 7. CREATE CLEAN ENSEMBL COLUMN
# ============================================================

df <- df %>%

  mutate(

    ENSEMBL = sub(
      "\\..*$",
      "",
      as.character(
        .data[[gene_col]]
      )
    )
  )


# ============================================================
# 8. ADD MOUSE GENE SYMBOLS
# ============================================================

df$SYMBOL <- suppressMessages(

  AnnotationDbi::mapIds(

    org.Mm.eg.db,

    keys = df$ENSEMBL,

    keytype = "ENSEMBL",

    column = "SYMBOL",

    multiVals = "first"
  )
)


df$SYMBOL <- as.character(
  df$SYMBOL
)


cat(
  "\nGenes with SYMBOL:",
  sum(
    !is.na(df$SYMBOL)
  ),
  "/",
  nrow(df),
  "\n"
)


# ============================================================
# 9. CREATE PRE-RANKED TABLE
#
# Ranking metric = DESeq2 Wald statistic
#
# stat > 0 = higher in HFpEF WT
# stat < 0 = higher in Control WT
# ============================================================

rank_df <- df %>%

  filter(
    !is.na(SYMBOL),
    SYMBOL != "",
    !is.na(stat)
  ) %>%

  mutate(
    rank_score = stat
  )


# ============================================================
# 10. REMOVE DUPLICATE SYMBOLS
# ============================================================

rank_df <- rank_df %>%

  arrange(
    desc(
      abs(rank_score)
    )
  ) %>%

  distinct(
    SYMBOL,
    .keep_all = TRUE
  )


cat(
  "\nUnique ranked gene symbols:",
  nrow(rank_df),
  "\n"
)


# ============================================================
# 11. CREATE GENE LIST
# ============================================================

gene_list <- rank_df$rank_score

names(gene_list) <- rank_df$SYMBOL


gene_list <- sort(
  gene_list,
  decreasing = TRUE
)


cat(
  "Maximum rank score:",
  max(gene_list),
  "\n"
)


cat(
  "Minimum rank score:",
  min(gene_list),
  "\n"
)


# ============================================================
# 12. HANDLE TIES
# ============================================================

n_ties <- sum(
  duplicated(gene_list)
)


cat(
  "Tied rank scores:",
  n_ties,
  "\n"
)


if (n_ties > 0) {

  tie_adjustment <- seq_along(
    gene_list
  ) * 1e-12


  gene_list <- gene_list +
    tie_adjustment


  gene_list <- sort(
    gene_list,
    decreasing = TRUE
  )
}


# ============================================================
# 13. EXPORT PRE-RANKED LIST
# ============================================================

rank_export <- tibble(

  SYMBOL = names(gene_list),

  rank_score = as.numeric(
    gene_list
  )
)


write_csv(
  rank_export,
  file.path(
    out_dir,
    "01_GSEA_preranked_gene_list.csv"
  )
)


write.table(
  rank_export,
  file = file.path(
    out_dir,
    "01_GSEA_preranked_gene_list.rnk"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)


# ============================================================
# 14. READ GLOBAL GO GMT
# ============================================================

gmt_all <- clusterProfiler::read.gmt(
  gmt_file
)


cat(
  "\nTotal GO gene sets in GMT:",
  length(
    unique(gmt_all$term)
  ),
  "\n"
)


cat(
  "Total TERM-GENE pairs:",
  nrow(gmt_all),
  "\n"
)


# ============================================================
# 15. SPLIT GO BEFORE GSEA
#
# IMPORTANT:
# BP / CC / MF are separated BEFORE multiple testing correction.
# Therefore each ontology gets its own BH/FDR correction.
# ============================================================

gmt_bp <- gmt_all %>%

  filter(
    str_detect(
      term,
      fixed("%GOBP%")
    )
  )


gmt_cc <- gmt_all %>%

  filter(
    str_detect(
      term,
      fixed("%GOCC%")
    )
  )


gmt_mf <- gmt_all %>%

  filter(
    str_detect(
      term,
      fixed("%GOMF%")
    )
  )


cat(
  "\nGO:BP gene sets:",
  length(
    unique(gmt_bp$term)
  ),
  "\n"
)


cat(
  "GO:CC gene sets:",
  length(
    unique(gmt_cc$term)
  ),
  "\n"
)


cat(
  "GO:MF gene sets:",
  length(
    unique(gmt_mf$term)
  ),
  "\n"
)


if (nrow(gmt_bp) == 0) {
  stop("No GOBP terms detected in GMT.")
}


if (nrow(gmt_cc) == 0) {
  stop("No GOCC terms detected in GMT.")
}


if (nrow(gmt_mf) == 0) {
  stop("No GOMF terms detected in GMT.")
}


# ============================================================
# 16. GSEA FUNCTION
# ============================================================

run_gsea <- function(gmt_subset) {

  clusterProfiler::GSEA(

    geneList = gene_list,

    TERM2GENE = gmt_subset,

    minGSSize = min_gs_size,

    maxGSSize = max_gs_size,

    pvalueCutoff = 1,

    pAdjustMethod = "BH",

    verbose = FALSE,

    seed = TRUE,

    by = "fgsea"
  )
}


# ============================================================
# 17. RUN GSEA BP / CC / MF
# ============================================================

set.seed(12345)

cat("\nRunning GO:BP GSEA...\n")

gsea_bp <- run_gsea(
  gmt_bp
)


set.seed(12345)

cat("Running GO:CC GSEA...\n")

gsea_cc <- run_gsea(
  gmt_cc
)


set.seed(12345)

cat("Running GO:MF GSEA...\n")

gsea_mf <- run_gsea(
  gmt_mf
)


# ============================================================
# 18. CONVERT RESULTS TO DATA FRAMES
# ============================================================

gsea_bp_df <- as.data.frame(
  gsea_bp
)


gsea_cc_df <- as.data.frame(
  gsea_cc
)


gsea_mf_df <- as.data.frame(
  gsea_mf
)


if (nrow(gsea_bp_df) == 0) {
  stop("GO:BP GSEA returned zero pathways.")
}


if (nrow(gsea_cc_df) == 0) {
  stop("GO:CC GSEA returned zero pathways.")
}


if (nrow(gsea_mf_df) == 0) {
  stop("GO:MF GSEA returned zero pathways.")
}


cat(
  "\nGO:BP pathways tested:",
  nrow(gsea_bp_df),
  "\n"
)


cat(
  "GO:CC pathways tested:",
  nrow(gsea_cc_df),
  "\n"
)


cat(
  "GO:MF pathways tested:",
  nrow(gsea_mf_df),
  "\n"
)


# ============================================================
# 19. ENSURE DESCRIPTION EXISTS
# ============================================================

if (!"Description" %in% names(gsea_bp_df)) {
  gsea_bp_df$Description <- gsea_bp_df$ID
}


if (!"Description" %in% names(gsea_cc_df)) {
  gsea_cc_df$Description <- gsea_cc_df$ID
}


if (!"Description" %in% names(gsea_mf_df)) {
  gsea_mf_df$Description <- gsea_mf_df$ID
}


# ============================================================
# 20. EXPORT ALL GSEA RESULTS
# ============================================================

write_csv(
  gsea_bp_df,
  file.path(
    out_dir,
    "02_GSEA_BP_all.csv"
  )
)


write_csv(
  gsea_cc_df,
  file.path(
    out_dir,
    "03_GSEA_CC_all.csv"
  )
)


write_csv(
  gsea_mf_df,
  file.path(
    out_dir,
    "04_GSEA_MF_all.csv"
  )
)


# ============================================================
# 21. SIGNIFICANT RESULTS
# ============================================================

gsea_bp_sig <- gsea_bp_df %>%

  filter(
    !is.na(p.adjust),
    p.adjust < gsea_fdr_cutoff
  ) %>%

  arrange(
    p.adjust
  )


gsea_cc_sig <- gsea_cc_df %>%

  filter(
    !is.na(p.adjust),
    p.adjust < gsea_fdr_cutoff
  ) %>%

  arrange(
    p.adjust
  )


gsea_mf_sig <- gsea_mf_df %>%

  filter(
    !is.na(p.adjust),
    p.adjust < gsea_fdr_cutoff
  ) %>%

  arrange(
    p.adjust
  )


write_csv(
  gsea_bp_sig,
  file.path(
    out_dir,
    "05_GSEA_BP_significant_FDR025.csv"
  )
)


write_csv(
  gsea_cc_sig,
  file.path(
    out_dir,
    "06_GSEA_CC_significant_FDR025.csv"
  )
)


write_csv(
  gsea_mf_sig,
  file.path(
    out_dir,
    "07_GSEA_MF_significant_FDR025.csv"
  )
)


cat(
  "\nSignificant GO:BP FDR < 0.25:",
  nrow(gsea_bp_sig),
  "\n"
)


cat(
  "Significant GO:CC FDR < 0.25:",
  nrow(gsea_cc_sig),
  "\n"
)


cat(
  "Significant GO:MF FDR < 0.25:",
  nrow(gsea_mf_sig),
  "\n"
)


# ============================================================
# 22. SPLIT UP / DOWN FUNCTION
# ============================================================

split_direction <- function(gsea_sig_df) {

  list(

    up = gsea_sig_df %>%
      filter(
        NES > 0
      ) %>%
      arrange(
        p.adjust
      ),

    down = gsea_sig_df %>%
      filter(
        NES < 0
      ) %>%
      arrange(
        p.adjust
      )
  )
}


bp_direction <- split_direction(
  gsea_bp_sig
)


cc_direction <- split_direction(
  gsea_cc_sig
)


mf_direction <- split_direction(
  gsea_mf_sig
)


bp_up <- bp_direction$up
bp_down <- bp_direction$down

cc_up <- cc_direction$up
cc_down <- cc_direction$down

mf_up <- mf_direction$up
mf_down <- mf_direction$down


cat(
  "\nGO:BP UP:",
  nrow(bp_up),
  " | DOWN:",
  nrow(bp_down),
  "\n"
)


cat(
  "GO:CC UP:",
  nrow(cc_up),
  " | DOWN:",
  nrow(cc_down),
  "\n"
)


cat(
  "GO:MF UP:",
  nrow(mf_up),
  " | DOWN:",
  nrow(mf_down),
  "\n"
)


# ============================================================
# 23. EXPORT UP / DOWN TABLES
# ============================================================

write_csv(
  bp_up,
  file.path(
    out_dir,
    "08_GSEA_BP_UP.csv"
  )
)


write_csv(
  bp_down,
  file.path(
    out_dir,
    "09_GSEA_BP_DOWN.csv"
  )
)


write_csv(
  cc_up,
  file.path(
    out_dir,
    "10_GSEA_CC_UP.csv"
  )
)


write_csv(
  cc_down,
  file.path(
    out_dir,
    "11_GSEA_CC_DOWN.csv"
  )
)


write_csv(
  mf_up,
  file.path(
    out_dir,
    "12_GSEA_MF_UP.csv"
  )
)


write_csv(
  mf_down,
  file.path(
    out_dir,
    "13_GSEA_MF_DOWN.csv"
  )
)


# ============================================================
# 24. GO DOTPLOT FUNCTION
#
# Selection:
#   - top positive NES terms by FDR
#   - top negative NES terms by FDR
#
# Display:
#   - positive and negative NES terms are grouped separately
#   - within each sign, pathways are ordered by gene set size
#   - largest gene set appears highest within its sign group
#   - point size = gene set size
#   - point colour = FDR (dark colours only)
# ============================================================

make_go_dotplot <- function(
  up_df,
  down_df,
  ontology,
  filename,
  highlight_lipid = FALSE,
  blue_category = c("none", "amino", "mito")
) {

  blue_category <- match.arg(blue_category)

  # ----------------------------------------------------------
  # DOUBLE SELECTION:
  # top pathways are selected independently in each direction
  # using the FDR-sorted UP and DOWN tables.
  # ----------------------------------------------------------

  top_up <- up_df %>%
    slice_head(n = n_top_up)

  top_down <- down_df %>%
    slice_head(n = n_top_down)

  top_df <- bind_rows(
    top_down,
    top_up
  )

  if (nrow(top_df) == 0) {
    message(
      "No significant ",
      ontology,
      " terms available for dotplot."
    )
    return(NULL)
  }

  top_df <- top_df %>%
    mutate(

      Description_clean = str_remove(
        Description,
        "%.*$"
      ),

      is_lipid = str_detect(
        Description_clean,
        regex(
          "fatty acid|lipid|lipoprotein|acyl|cholesterol|triglycer|phospholipid|PPAR",
          ignore_case = TRUE
        )
      ),

      is_amino = str_detect(
        Description_clean,
        regex(
          "amino acid|branched.chain|valine|leucine|isoleucine|alanine|aspartate|glutamate|2-oxocarboxylic",
          ignore_case = TRUE
        )
      ),

      is_mito = str_detect(
        Description_clean,
        regex(
          "mitochond|oxidative phosphorylation|respiratory chain",
          ignore_case = TRUE
        )
      ),

      is_blue = case_when(
        blue_category == "amino" ~ is_amino,
        blue_category == "mito"  ~ is_mito,
        TRUE ~ FALSE
      ),

      Description_label = case_when(

        highlight_lipid & is_lipid ~ paste0(
          "<span style='color:#D55E00'>",
          Description_clean,
          "</span>"
        ),

        is_blue ~ paste0(
          "<span style='color:#0072B2'>",
          Description_clean,
          "</span>"
        ),

        TRUE ~ Description_clean
      ),

      direction = if_else(
        NES < 0,
        "Negative NES",
        "Positive NES"
      )
    )

  # ----------------------------------------------------------
  # Y-axis ordering
  #
  # Negative pathways form one block and positive pathways
  # form another block. Within each block, gene set size is
  # used for ordering. Because the last factor levels appear
  # at the top in ggplot, the positive block is displayed on
  # top and the largest set is highest within each block.
  # ----------------------------------------------------------

  top_df <- top_df %>%
    mutate(
      direction = factor(
        direction,
        levels = c(
          "Negative NES",
          "Positive NES"
        )
      )
    ) %>%
    arrange(
      direction,
      setSize
    ) %>%
    mutate(
      Description_plot = factor(
        Description_label,
        levels = Description_label
      )
    )

  n_negative <- sum(
    top_df$NES < 0
  )

  nes_limit <- ceiling(
    max(
      abs(top_df$NES),
      na.rm = TRUE
    ) * 10
  ) / 10

  p <- ggplot(
    top_df,
    aes(
      x = NES,
      y = Description_plot
    )
  ) +

    geom_point(
      aes(
        size = setSize,
        colour = p.adjust
      )
    ) +

    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5
    ) +

    {
      if (
        n_negative > 0 &&
        n_negative < nrow(top_df)
      ) {
        geom_hline(
          yintercept = n_negative + 0.5,
          colour = "grey55",
          linewidth = 0.3
        )
      }
    } +

    scale_x_continuous(
      limits = c(
        -nes_limit,
        nes_limit
      ),
      breaks = scales::breaks_pretty(
        n = 7
      ),
      labels = scales::label_number(
        accuracy = 0.1
      ),
      expand = expansion(
        mult = c(
          0.05,
          0.05
        )
      )
    ) +

    scale_size_continuous(
      range = c(
        2.5,
        8
      ),
      breaks = scales::breaks_pretty(
        n = 4
      )
    ) +

    scale_colour_gradientn(
      colours = c(
        "#081D58",
        "#253494",
        "#54278F",
        "#8C2D75",
        "#B2182B"
      ),
      trans = "log10"
    ) +

    theme_classic(
      base_size = 11
    ) +

    theme(
      axis.text.y = ggtext::element_markdown(
        size = 9
      ),
      axis.text.x = element_text(
        size = 10
      ),
      axis.title.x = element_text(
        size = 11
      ),
      plot.title = element_text(
        hjust = 0.5
      ),
      plot.margin = margin(
        t = 10,
        r = 20,
        b = 10,
        l = 10
      )
    ) +

    labs(
      title = paste0(
        "GSEA ",
        ontology
      ),
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      size = "Gene set size",
      colour = "FDR"
    )

  ggsave(
    file.path(
      out_dir,
      paste0(
        filename,
        ".png"
      )
    ),
    p,
    width = 10,
    height = 7,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    file.path(
      out_dir,
      paste0(
        filename,
        ".pdf"
      )
    ),
    p,
    width = 10,
    height = 7,
    bg = "white"
  )

  return(p)
}


# ============================================================
# 25A. PANEL C -- GO:BP GLOBAL
#
# Lipid-related labels      = orange
# Amino-acid-related labels = blue
# All other labels          = black
# ============================================================

p_bp <- make_go_dotplot(
  bp_up,
  bp_down,
  "GO:BP",
  "14_GSEA_BP_top_pathways_dotplot",
  highlight_lipid = TRUE,
  blue_category = "amino"
)


# ============================================================
# 25B. GENERAL GO:CC
#
# Mitochondrial labels = blue
# Exploratory / supplementary global CC view.
# ============================================================

p_cc <- make_go_dotplot(
  cc_up,
  cc_down,
  "GO:CC",
  "15_GSEA_CC_top_pathways_dotplot",
  highlight_lipid = FALSE,
  blue_category = "mito"
)


# ============================================================
# 25C. GO:MF
# ============================================================

p_mf <- make_go_dotplot(
  mf_up,
  mf_down,
  "GO:MF",
  "16_GSEA_MF_top_pathways_dotplot",
  highlight_lipid = FALSE,
  blue_category = "none"
)

# ============================================================
# 25D. PANEL E -- TOP NEGATIVELY ENRICHED GO:CC
#
# IMPORTANT:
# Selection is NOT based on mitochondrial annotation.
# 1) significant GO:CC only
# 2) NES < 0
# 3) strongest negative NES
# Mitochondrial terms are highlighted only AFTER selection.
# ============================================================

cc_negative_top <- gsea_cc_sig %>%
  filter(
    NES < 0
  ) %>%
  arrange(
    NES
  ) %>%
  slice_head(
    n = n_top_cc_negative
  ) %>%
  mutate(

    Description_clean = str_remove(
      Description,
      "%.*$"
    ),

    is_mito = str_detect(
      Description_clean,
      regex(
        "mitochond",
        ignore_case = TRUE
      )
    ),

    Description_label = if_else(
      is_mito,
      paste0(
        "<span style='color:#0072B2'>",
        Description_clean,
        "</span>"
      ),
      Description_clean
    )
  )


write_csv(
  cc_negative_top,
  file.path(
    out_dir,
    "17_GSEA_CC_top_negative.csv"
  )
)


if (nrow(cc_negative_top) > 0) {

  cc_negative_top <- cc_negative_top %>%
    arrange(
      setSize
    ) %>%
    mutate(
      Description_plot = factor(
        Description_label,
        levels = Description_label
      )
    )


  p_cc_negative <- ggplot(
    cc_negative_top,
    aes(
      x = NES,
      y = Description_plot
    )
  ) +

    geom_point(
      aes(
        size = setSize,
        colour = p.adjust
      )
    ) +

    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5
    ) +

    scale_x_continuous(
      breaks = scales::breaks_pretty(
        n = 6
      ),
      labels = scales::label_number(
        accuracy = 0.1
      )
    ) +

    scale_size_continuous(
      range = c(
        2.5,
        8
      ),
      breaks = scales::breaks_pretty(
        n = 4
      )
    ) +

    scale_colour_gradientn(
      colours = c(
        "#081D58",
        "#253494",
        "#54278F",
        "#8C2D75",
        "#B2182B"
      ),
      trans = "log10"
    ) +

    theme_classic(
      base_size = 11
    ) +

    theme(

      axis.text.y = ggtext::element_markdown(
        size = 9
      ),

      axis.text.x = element_text(
        size = 10
      ),

      axis.title.x = element_text(
        size = 11
      ),

      plot.title = element_text(
        hjust = 0.5
      )
    ) +

    labs(
      title = "Negatively enriched GO:CC",
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      size = "Gene set size",
      colour = "FDR"
    )


  ggsave(
    file.path(
      out_dir,
      "17_GSEA_CC_top_negative_dotplot.png"
    ),
    p_cc_negative,
    width = 10,
    height = 6,
    dpi = 600,
    bg = "white"
  )


  ggsave(
    file.path(
      out_dir,
      "17_GSEA_CC_top_negative_dotplot.pdf"
    ),
    p_cc_negative,
    width = 10,
    height = 6,
    bg = "white"
  )
}


# ============================================================
# 26. MITOCHONDRIAL GO:CC SEARCH
#
# Search ALL CC results, not only significant ones.
# ============================================================

mito_cc <- gsea_cc_df %>%

  filter(

    str_detect(
      Description,
      regex(
        "mitochond",
        ignore_case = TRUE
      )
    )
  ) %>%

  arrange(
    p.adjust
  )


write_csv(
  mito_cc,
  file.path(
    out_dir,
    "SUPP_GSEA_CC_mitochondrial_all.csv"
  )
)


mito_cc_sig <- mito_cc %>%

  filter(
    !is.na(p.adjust),
    p.adjust < gsea_fdr_cutoff
  )


write_csv(
  mito_cc_sig,
  file.path(
    out_dir,
    "SUPP_GSEA_CC_mitochondrial_FDR025.csv"
  )
)


cat(
  "\nMitochondrial GO:CC terms detected:",
  nrow(mito_cc),
  "\n"
)


cat(
  "Significant mitochondrial GO:CC terms:",
  nrow(mito_cc_sig),
  "\n"
)


if (nrow(mito_cc) > 0) {

  cat(
    "\nTop mitochondrial GO:CC terms:\n"
  )


  print(

    mito_cc %>%

      dplyr::select(
        Description,
        setSize,
        NES,
        pvalue,
        p.adjust
      ) %>%

      slice_head(
        n = 20
      )
  )
}


# ============================================================
# 27. MITOCHONDRIAL CC DOTPLOT
# ============================================================

if (nrow(mito_cc_sig) > 0) {

  mito_cc_top <- mito_cc_sig %>%

    arrange(
      p.adjust,
      desc(
        abs(NES)
      )
    ) %>%

    slice_head(
      n = 20
    ) %>%

    mutate(

      Description_clean = str_remove(
        Description,
        "%.*$"
      ),

      Description_plot = factor(
        Description_clean,
        levels = Description_clean[
          order(setSize)
        ]
      )
    )


  p_mito_cc <- ggplot(

    mito_cc_top,

    aes(
      x = NES,
      y = Description_plot
    )
  ) +

    geom_point(

      aes(
        size = setSize,
        colour = p.adjust
      )
    ) +

    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5
    ) +

    scale_x_continuous(
      breaks = scales::breaks_pretty(
        n = 6
      ),
      labels = scales::label_number(
        accuracy = 0.1
      )
    ) +

    scale_size_continuous(
      range = c(
        2.5,
        8
      ),
      breaks = scales::breaks_pretty(
        n = 4
      )
    ) +

    scale_colour_gradientn(
      colours = c(
        "#081D58",
        "#253494",
        "#54278F",
        "#8C2D75",
        "#B2182B"
      ),
      trans = "log10"
    ) +

    theme_classic(
      base_size = 11
    ) +

    labs(
      title = "Mitochondrial GO:CC",
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      size = "Gene set size",
      colour = "FDR"
    )


  ggsave(
    file.path(
      out_dir,
      "SUPP_GSEA_CC_mitochondrial_dotplot.png"
    ),
    p_mito_cc,
    width = 9,
    height = 6,
    dpi = 600,
    bg = "white"
  )


  ggsave(
    file.path(
      out_dir,
      "SUPP_GSEA_CC_mitochondrial_dotplot.pdf"
    ),
    p_mito_cc,
    width = 9,
    height = 6,
    bg = "white"
  )
}


# ============================================================
# 28. BP METABOLIC PATHWAY SEARCH
#
# Biological metabolic processes belong mainly to GO:BP.
# ============================================================

metabolic_pattern <- paste(

  c(
    "fatty acid",
    "lipid",
    "PPAR",
    "peroxisome",
    "glucose",
    "glycol",
    "pyruvate",
    "amino acid",
    "branched.chain",
    "leucine",
    "isoleucine",
    "valine",
    "ketone",
    "ketogenesis",
    "oxidative phosphorylation",
    "mitochondrial",
    "acetyl.CoA",
    "TCA",
    "citric acid",
    "tricarboxylic"
  ),

  collapse = "|"
)


gsea_metabolic <- gsea_bp_df %>%

  filter(

    str_detect(

      Description,

      regex(
        metabolic_pattern,
        ignore_case = TRUE
      )
    )
  ) %>%

  arrange(
    p.adjust
  )


write_csv(
  gsea_metabolic,
  file.path(
    out_dir,
    "20_GSEA_BP_metabolic_all.csv"
  )
)


gsea_metabolic_sig <- gsea_metabolic %>%

  filter(
    !is.na(p.adjust),
    p.adjust < gsea_fdr_cutoff
  )


write_csv(
  gsea_metabolic_sig,
  file.path(
    out_dir,
    "21_GSEA_BP_metabolic_FDR025.csv"
  )
)


cat(
  "\nSignificant metabolic GO:BP pathways:",
  nrow(gsea_metabolic_sig),
  "\n"
)


# ============================================================
# 29. METABOLIC BP DOTPLOT
#
# Supplementary targeted view.
# Positive and negative NES terms are grouped separately.
# ============================================================

metabolic_top <- gsea_metabolic_sig %>%
  arrange(
    p.adjust,
    desc(
      abs(NES)
    )
  ) %>%
  slice_head(
    n = 25
  )


if (nrow(metabolic_top) > 0) {

  metabolic_top <- metabolic_top %>%
    mutate(
      Description_clean = str_remove(
        Description,
        "%.*$"
      ),
      direction = if_else(
        NES < 0,
        "Negative NES",
        "Positive NES"
      ),
      direction = factor(
        direction,
        levels = c(
          "Negative NES",
          "Positive NES"
        )
      )
    ) %>%
    arrange(
      direction,
      setSize
    ) %>%
    mutate(
      Description_plot = factor(
        Description_clean,
        levels = Description_clean
      )
    )

  n_negative_metabolic <- sum(
    metabolic_top$NES < 0
  )

  nes_limit_metabolic <- ceiling(
    max(
      abs(metabolic_top$NES),
      na.rm = TRUE
    ) * 10
  ) / 10

  p_metabolic <- ggplot(
    metabolic_top,
    aes(
      x = NES,
      y = Description_plot
    )
  ) +

    geom_point(
      aes(
        size = setSize,
        colour = p.adjust
      )
    ) +

    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +

    {
      if (
        n_negative_metabolic > 0 &&
        n_negative_metabolic < nrow(metabolic_top)
      ) {
        geom_hline(
          yintercept = n_negative_metabolic + 0.5,
          colour = "grey55",
          linewidth = 0.3
        )
      }
    } +

    scale_x_continuous(
      limits = c(
        -nes_limit_metabolic,
        nes_limit_metabolic
      ),
      breaks = scales::breaks_pretty(
        n = 7
      ),
      labels = scales::label_number(
        accuracy = 0.1
      )
    ) +

    scale_size_continuous(
      range = c(
        2.5,
        8
      ),
      breaks = scales::breaks_pretty(
        n = 4
      )
    ) +

    scale_colour_gradientn(
      colours = c(
        "#081D58",
        "#253494",
        "#54278F",
        "#8C2D75",
        "#B2182B"
      ),
      trans = "log10"
    ) +

    theme_classic(
      base_size = 11
    ) +

    labs(
      title = "Metabolic GO:BP",
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      size = "Gene set size",
      colour = "FDR"
    )

  ggsave(
    file.path(
      out_dir,
      "22_GSEA_BP_metabolic_dotplot.png"
    ),
    p_metabolic,
    width = 8,
    height = 7,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    file.path(
      out_dir,
      "22_GSEA_BP_metabolic_dotplot.pdf"
    ),
    p_metabolic,
    width = 8,
    height = 7,
    bg = "white"
  )
}

# ============================================================
# 30. FATTY ACID METABOLIC PROCESS
#
# GO:BP ONLY
# ============================================================

fatty_hits <- gsea_bp_df %>%

  filter(

    str_detect(

      Description,

      regex(
        "^FATTY ACID METABOLIC PROCESS%GOBP%",
        ignore_case = TRUE
      )
    )
  ) %>%

  arrange(
    p.adjust
  )


write_csv(
  fatty_hits,
  file.path(
    out_dir,
    "23_GSEA_BP_fatty_acid_hits.csv"
  )
)


cat(
  "\nExact fatty acid metabolic process hits:",
  nrow(fatty_hits),
  "\n"
)


if (nrow(fatty_hits) > 0) {

  print(

    fatty_hits %>%

      dplyr::select(
        ID,
        Description,
        setSize,
        enrichmentScore,
        NES,
        pvalue,
        p.adjust
      )
  )
}


# ============================================================
# 31. FATTY ACID GSEA CURVE
# ============================================================

if (nrow(fatty_hits) > 0) {

  fatty_id <- fatty_hits$ID[1]


  fatty_name <- str_remove(
    fatty_hits$Description[1],
    "%.*$"
  )


  fatty_name <- str_to_sentence(
    tolower(
      fatty_name
    )
  )


  fatty_nes <- round(
    fatty_hits$NES[1],
    2
  )


  fatty_fdr <- formatC(
    fatty_hits$p.adjust[1],
    format = "f",
    digits = 6
  )


  p_fatty <- enrichplot::gseaplot2(

    gsea_bp,

    geneSetID = fatty_id,

    title = paste0(
      "[GO:BP] ",
      fatty_name
    ),

    pvalue_table = FALSE,

    ES_geom = "line",

    subplots = 1:3,

    rel_heights = c(
      1.6,
      0.45,
      0.95
    ),

    base_size = 11
  )


  clean_panel <- theme(

    panel.grid.major = element_blank(),

    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.6
    ),

    axis.line = element_blank(),

    axis.text = element_text(
      size = 9
    ),

    axis.title = element_text(
      size = 10
    ),

    plot.background = element_rect(
      fill = "white",
      colour = NA
    )
  )


  p_fatty[[1]] <- p_fatty[[1]] +

    labs(
      y = "Enrichment Score"
    ) +

    clean_panel +

    theme(

      plot.title = element_text(
        hjust = 0.5,
        size = 12,
        face = "plain"
      ),

      plot.margin = margin(
        t = 5,
        r = 5,
        b = 0,
        l = 5
      )
    ) +

    annotate(

      "text",

      x = Inf,

      y = Inf,

      label = paste0(
        "NES = ",
        fatty_nes,
        "\nFDR = ",
        fatty_fdr
      ),

      hjust = 1.08,

      vjust = 1.15,

      size = 3.8
    )


  p_fatty[[2]] <- p_fatty[[2]] +

    clean_panel +

    theme(

      axis.title = element_blank(),

      axis.text.y = element_blank(),

      axis.ticks.y = element_blank(),

      plot.margin = margin(
        t = 0,
        r = 5,
        b = 0,
        l = 5
      )
    )


  p_fatty[[3]] <- p_fatty[[3]] +

    labs(
      x = "Rank in Ordered Dataset",
      y = "Ranked list metric"
    ) +

    clean_panel +

    theme(

      plot.margin = margin(
        t = 0,
        r = 5,
        b = 5,
        l = 5
      )
    )


  ggsave(
    file.path(
      out_dir,
      "24_GSEA_BP_fatty_acid_curve.png"
    ),
    plot = p_fatty,
    width = 5.2,
    height = 4.2,
    dpi = 600,
    bg = "white"
  )


  ggsave(
    file.path(
      out_dir,
      "24_GSEA_BP_fatty_acid_curve.pdf"
    ),
    plot = p_fatty,
    width = 5.2,
    height = 4.2,
    bg = "white"
  )
}


# ============================================================
# 32. CORE / LEADING-EDGE GENES
#
# BP
# ============================================================

if ("core_enrichment" %in% names(gsea_bp_df)) {

  core_bp <- gsea_bp_df %>%

    filter(
      !is.na(core_enrichment),
      core_enrichment != ""
    ) %>%

    dplyr::select(
      ID,
      Description,
      NES,
      pvalue,
      p.adjust,
      core_enrichment
    ) %>%

    separate_rows(
      core_enrichment,
      sep = "/"
    ) %>%

    rename(
      SYMBOL = core_enrichment
    )


  write_csv(
    core_bp,
    file.path(
      out_dir,
      "25_GSEA_BP_core_enrichment_genes.csv"
    )
  )


  core_metabolic <- core_bp %>%

    filter(

      str_detect(

        Description,

        regex(
          metabolic_pattern,
          ignore_case = TRUE
        )
      )
    )


  write_csv(
    core_metabolic,
    file.path(
      out_dir,
      "26_GSEA_BP_metabolic_core_genes.csv"
    )
  )
}

# ============================================================
# 33. PANEL F -- KEGG ORA
#
# Over-representation analysis on DEGs
# DEG criterion: padj < 0.05
#
# No log2FC cutoff.
# Global KEGG analysis: no metabolic pre-selection.
# ============================================================

# ------------------------------------------------------------
# Check DESeq2 adjusted p-values
# ------------------------------------------------------------

if (!"padj" %in% names(df)) {
  stop("Column 'padj' is required for KEGG ORA.")
}


# ------------------------------------------------------------
# ENSEMBL -> ENTREZ
# ------------------------------------------------------------

df$ENTREZID <- suppressMessages(

  AnnotationDbi::mapIds(

    org.Mm.eg.db,

    keys = df$ENSEMBL,

    keytype = "ENSEMBL",

    column = "ENTREZID",

    multiVals = "first"
  )
)


df$ENTREZID <- as.character(
  df$ENTREZID
)


# ------------------------------------------------------------
# DEG LIST
#
# padj < 0.05
# No fold-change cutoff
# ------------------------------------------------------------

kegg_deg <- df %>%

  filter(
    !is.na(padj),
    padj < 0.05,
    !is.na(ENTREZID),
    ENTREZID != ""
  ) %>%

  arrange(
    padj
  ) %>%

  distinct(
    ENTREZID,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# BACKGROUND
#
# Genes eligible for the DESeq2 significance test
# and successfully mapped to ENTREZ.
# ------------------------------------------------------------

kegg_background <- df %>%

  filter(
    !is.na(padj),
    !is.na(ENTREZID),
    ENTREZID != ""
  ) %>%

  distinct(
    ENTREZID
  ) %>%

  pull(
    ENTREZID
  )


cat(
  "\n==============================\n",
  "KEGG ORA\n",
  "==============================\n",
  sep = ""
)


cat(
  "DEGs used for KEGG ORA:",
  nrow(kegg_deg),
  "\n"
)


cat(
  "Background genes:",
  length(kegg_background),
  "\n"
)


# ------------------------------------------------------------
# RUN KEGG ORA
# ------------------------------------------------------------

kegg_ora <- clusterProfiler::enrichKEGG(

  gene = kegg_deg$ENTREZID,

  universe = kegg_background,

  organism = "mmu",

  keyType = "ncbi-geneid",

  pvalueCutoff = 1,

  pAdjustMethod = "BH",

  qvalueCutoff = 1
)


kegg_ora_df <- as.data.frame(
  kegg_ora
)


if (nrow(kegg_ora_df) == 0) {

  warning(
    "KEGG ORA returned zero pathways."
  )

} else {

  kegg_ora_df <- kegg_ora_df %>%

    arrange(
      p.adjust
    )
}


# ------------------------------------------------------------
# EXPORT ALL RESULTS
# ------------------------------------------------------------

write_csv(

  kegg_ora_df,

  file.path(
    out_dir,
    "27_KEGG_ORA_all.csv"
  )
)


# ------------------------------------------------------------
# Significant results
#
# Here use conventional ORA threshold FDR < 0.05
# NOT GSEA FDR < 0.25
# ------------------------------------------------------------

kegg_ora_sig <- kegg_ora_df %>%

  filter(
    !is.na(p.adjust),
    p.adjust < 0.05
  ) %>%

  arrange(
    p.adjust
  )


write_csv(

  kegg_ora_sig,

  file.path(
    out_dir,
    "28_KEGG_ORA_significant_FDR005.csv"
  )
)


cat(
  "KEGG pathways tested:",
  nrow(kegg_ora_df),
  "\n"
)


cat(
  "Significant KEGG pathways FDR < 0.05:",
  nrow(kegg_ora_sig),
  "\n"
)


# ============================================================
# 34. PANEL F -- GLOBAL KEGG ORA DOTPLOT
#
# IMPORTANT:
# Top pathways selected globally by adjusted p-value.
# No biological / metabolic filtering.
#
# x      = GeneRatio
# size   = Count
# colour = FDR
# ============================================================

n_top_kegg_ora <- 20


if (nrow(kegg_ora_df) > 0) {

  kegg_ora_top <- kegg_ora_df %>%

    filter(
      !is.na(p.adjust)
    ) %>%

    arrange(
      p.adjust
    ) %>%

    slice_head(
      n = n_top_kegg_ora
    ) %>%

    mutate(

      GeneRatio_numeric =

        as.numeric(
          sub(
            "/.*",
            "",
            GeneRatio
          )
        ) /

        as.numeric(
          sub(
            ".*/",
            "",
            GeneRatio
          )
        )
    )


  # ----------------------------------------------------------
  # Display order:
  # largest Count at top
  # ----------------------------------------------------------

  kegg_ora_top <- kegg_ora_top %>%

    arrange(
      Count
    ) %>%

    mutate(

      Description_plot = factor(
        Description,
        levels = Description
      )
    )


  p_kegg_ora <- ggplot(

    kegg_ora_top,

    aes(
      x = GeneRatio_numeric,
      y = Description_plot
    )
  ) +

    geom_point(

      aes(
        size = Count,
        colour = p.adjust
      )
    ) +

    scale_size_continuous(
      range = c(
        2.5,
        8
      ),
      breaks = scales::breaks_pretty(
        n = 4
      )
    ) +

    scale_colour_gradientn(
      colours = c(
        "#081D58",
        "#253494",
        "#54278F",
        "#8C2D75",
        "#B2182B"
      ),
      trans = "log10"
    ) +

    scale_x_continuous(

      labels = scales::label_number(
        accuracy = 0.01
      )
    ) +

    theme_classic(
      base_size = 11
    ) +

    theme(

      axis.text.y = element_text(
        size = 9
      ),

      axis.text.x = element_text(
        size = 10
      ),

      axis.title.x = element_text(
        size = 11
      ),

      plot.title = element_text(
        hjust = 0.5
      ),

      plot.margin = margin(
        t = 10,
        r = 20,
        b = 10,
        l = 10
      )
    ) +

    labs(

      title = "KEGG pathway enrichment",

      x = "Gene ratio",

      y = NULL,

      size = "Gene count",

      colour = "FDR"
    )


  ggsave(

    file.path(
      out_dir,
      "29_KEGG_ORA_top_pathways_dotplot.png"
    ),

    p_kegg_ora,

    width = 10,

    height = 7,

    dpi = 600,

    bg = "white"
  )


  ggsave(

    file.path(
      out_dir,
      "29_KEGG_ORA_top_pathways_dotplot.pdf"
    ),

    p_kegg_ora,

    width = 10,

    height = 7,

    bg = "white"
  )
}


# ============================================================
# 35. PEROXISOME -- REPORT ONLY
#
# This does NOT affect pathway selection in the dotplot.
# It simply reports the Peroxisome statistics separately.
# ============================================================

kegg_ora_peroxisome <- kegg_ora_df %>%

  filter(
    ID == "mmu04146" |
      str_detect(
        Description,
        regex(
          "^Peroxisome$",
          ignore_case = TRUE
        )
      )
  )


write_csv(

  kegg_ora_peroxisome,

  file.path(
    out_dir,
    "30_KEGG_ORA_peroxisome.csv"
  )
)


if (nrow(kegg_ora_peroxisome) > 0) {

  cat(
    "\nKEGG ORA -- PEROXISOME\n"
  )

  print(

    as_tibble(
      kegg_ora_peroxisome %>%

        dplyr::select(
          ID,
          Description,
          GeneRatio,
          BgRatio,
          Count,
          pvalue,
          p.adjust
        )
    )
  )
}

# ============================================================
# 33. KEGG GSEA
#
# Same DESeq2 Wald-statistic ranking as GO GSEA.
# KEGG requires ENTREZ IDs.
# ============================================================

rank_kegg <- rank_df %>%
  dplyr::select(
    SYMBOL,
    rank_score
  )


rank_kegg$ENTREZID <- suppressMessages(
  AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = rank_kegg$SYMBOL,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )
)


rank_kegg <- rank_kegg %>%
  filter(
    !is.na(ENTREZID),
    ENTREZID != "",
    !is.na(rank_score)
  ) %>%
  arrange(
    desc(
      abs(rank_score)
    )
  ) %>%
  distinct(
    ENTREZID,
    .keep_all = TRUE
  )


gene_list_kegg <- rank_kegg$rank_score
names(gene_list_kegg) <- rank_kegg$ENTREZID


gene_list_kegg <- sort(
  gene_list_kegg,
  decreasing = TRUE
)


# Handle possible ties after SYMBOL -> ENTREZ conversion.
n_ties_kegg <- sum(
  duplicated(gene_list_kegg)
)


if (n_ties_kegg > 0) {

  tie_adjustment_kegg <- seq_along(
    gene_list_kegg
  ) * 1e-12

  gene_list_kegg <- gene_list_kegg +
    tie_adjustment_kegg

  gene_list_kegg <- sort(
    gene_list_kegg,
    decreasing = TRUE
  )
}


cat(
  "\nKEGG ranked ENTREZ genes:",
  length(gene_list_kegg),
  "\n"
)


set.seed(12345)

cat(
  "Running KEGG GSEA...\n"
)


gsea_kegg <- clusterProfiler::gseKEGG(
  geneList = gene_list_kegg,
  organism = "mmu",
  keyType = "ncbi-geneid",
  minGSSize = min_gs_size,
  maxGSSize = max_gs_size,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  seed = TRUE,
  by = "fgsea"
)


gsea_kegg_df <- as.data.frame(
  gsea_kegg
)


if (!"Description" %in% names(gsea_kegg_df) & nrow(gsea_kegg_df) > 0) {
  gsea_kegg_df$Description <- gsea_kegg_df$ID
}


write_csv(
  gsea_kegg_df,
  file.path(
    out_dir,
    "SUPP_GSEA_KEGG_all.csv"
  )
)


if (nrow(gsea_kegg_df) > 0) {

  gsea_kegg_sig <- gsea_kegg_df %>%
    filter(
      !is.na(p.adjust),
      p.adjust < gsea_fdr_cutoff
    ) %>%
    arrange(
      p.adjust
    )

} else {

  gsea_kegg_sig <- gsea_kegg_df
}


write_csv(
  gsea_kegg_sig,
  file.path(
    out_dir,
    "SUPP_GSEA_KEGG_significant_FDR025.csv"
  )
)


kegg_up <- gsea_kegg_sig %>%
  filter(
    NES > 0
  ) %>%
  arrange(
    p.adjust
  )


kegg_down <- gsea_kegg_sig %>%
  filter(
    NES < 0
  ) %>%
  arrange(
    p.adjust
  )


cat(
  "KEGG pathways tested:",
  nrow(gsea_kegg_df),
  "\n"
)

cat(
  "Significant KEGG FDR < 0.25:",
  nrow(gsea_kegg_sig),
  " | UP:",
  nrow(kegg_up),
  " | DOWN:",
  nrow(kegg_down),
  "\n"
)


# ============================================================
# 34. SUPPLEMENTARY KEGG GSEA DOTPLOT
#
# Top significant positive + negative pathways by FDR.
# Positive and negative NES terms are grouped separately.
# ============================================================

kegg_top <- bind_rows(
  kegg_down %>%
    slice_head(
      n = n_top_kegg_down
    ),
  kegg_up %>%
    slice_head(
      n = n_top_kegg_up
    )
)


if (nrow(kegg_top) > 0) {

  kegg_top <- kegg_top %>%
    mutate(
      direction = if_else(
        NES < 0,
        "Negative NES",
        "Positive NES"
      ),
      direction = factor(
        direction,
        levels = c(
          "Negative NES",
          "Positive NES"
        )
      )
    ) %>%
    arrange(
      direction,
      setSize
    ) %>%
    mutate(
      Description_plot = factor(
        Description,
        levels = Description
      )
    )

  n_negative_kegg <- sum(
    kegg_top$NES < 0
  )

  nes_limit_kegg <- ceiling(
    max(
      abs(kegg_top$NES),
      na.rm = TRUE
    ) * 10
  ) / 10

  p_kegg <- ggplot(
    kegg_top,
    aes(
      x = NES,
      y = Description_plot
    )
  ) +

    geom_point(
      aes(
        size = setSize,
        colour = p.adjust
      )
    ) +

    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5
    ) +

    {
      if (
        n_negative_kegg > 0 &&
        n_negative_kegg < nrow(kegg_top)
      ) {
        geom_hline(
          yintercept = n_negative_kegg + 0.5,
          colour = "grey55",
          linewidth = 0.3
        )
      }
    } +

    scale_x_continuous(
      limits = c(
        -nes_limit_kegg,
        nes_limit_kegg
      ),
      breaks = scales::breaks_pretty(
        n = 7
      ),
      labels = scales::label_number(
        accuracy = 0.1
      ),
      expand = expansion(
        mult = c(
          0.05,
          0.05
        )
      )
    ) +

    scale_size_continuous(
      range = c(
        2.5,
        8
      ),
      breaks = scales::breaks_pretty(
        n = 4
      )
    ) +

    scale_colour_gradientn(
      colours = c(
        "#081D58",
        "#253494",
        "#54278F",
        "#8C2D75",
        "#B2182B"
      ),
      trans = "log10"
    ) +

    theme_classic(
      base_size = 11
    ) +

    theme(
      axis.text.y = element_text(
        size = 9
      ),
      axis.text.x = element_text(
        size = 10
      ),
      axis.title.x = element_text(
        size = 11
      ),
      plot.title = element_text(
        hjust = 0.5
      )
    ) +

    labs(
      title = "GSEA KEGG",
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      size = "Gene set size",
      colour = "FDR"
    )

  ggsave(
    file.path(
      out_dir,
      "SUPP_GSEA_KEGG_top_pathways_dotplot.png"
    ),
    p_kegg,
    width = 10,
    height = 7,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    file.path(
      out_dir,
      "SUPP_GSEA_KEGG_top_pathways_dotplot.pdf"
    ),
    p_kegg,
    width = 10,
    height = 7,
    bg = "white"
  )
}

# ============================================================
# 35. KEGG PEROXISOME PATHWAY
#
# Mouse KEGG Peroxisome = mmu04146.
# We test what KEGG actually returned; significance is not assumed.
# ============================================================

if (nrow(gsea_kegg_df) > 0) {

  peroxisome_hits <- gsea_kegg_df %>%
    filter(
      ID == "mmu04146" |
        str_detect(
          Description,
          regex(
            "^Peroxisome$",
            ignore_case = TRUE
          )
        )
    )

} else {

  peroxisome_hits <- gsea_kegg_df
}


write_csv(
  peroxisome_hits,
  file.path(
    out_dir,
    "SUPP_GSEA_KEGG_peroxisome.csv"
  )
)


if (nrow(peroxisome_hits) > 0) {

  cat(
    "\nKEGG PEROXISOME\n"
  )

  print(
    peroxisome_hits %>%
      dplyr::select(
        ID,
        Description,
        setSize,
        NES,
        pvalue,
        p.adjust
      )
  )


  perox_id <- peroxisome_hits$ID[1]

  perox_nes <- round(
    peroxisome_hits$NES[1],
    2
  )

  perox_fdr <- formatC(
    peroxisome_hits$p.adjust[1],
    format = "f",
    digits = 6
  )


  p_perox <- enrichplot::gseaplot2(
    gsea_kegg,
    geneSetID = perox_id,
    title = "[KEGG] Peroxisome",
    pvalue_table = FALSE,
    ES_geom = "line",
    subplots = 1:3,
    rel_heights = c(
      1.6,
      0.45,
      0.95
    ),
    base_size = 11
  )


  clean_panel_kegg <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.6
    ),
    axis.line = element_blank(),
    axis.text = element_text(
      size = 9
    ),
    axis.title = element_text(
      size = 10
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    )
  )


  p_perox[[1]] <- p_perox[[1]] +
    labs(
      y = "Enrichment Score"
    ) +
    clean_panel_kegg +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 12,
        face = "plain"
      ),
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 0,
        l = 5
      )
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste0(
        "NES = ",
        perox_nes,
        "\nFDR = ",
        perox_fdr
      ),
      hjust = 1.08,
      vjust = 1.15,
      size = 3.8
    )


  p_perox[[2]] <- p_perox[[2]] +
    clean_panel_kegg +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(
        t = 0,
        r = 5,
        b = 0,
        l = 5
      )
    )


  p_perox[[3]] <- p_perox[[3]] +
    labs(
      x = "Rank in Ordered Dataset",
      y = "Ranked list metric"
    ) +
    clean_panel_kegg +
    theme(
      plot.margin = margin(
        t = 0,
        r = 5,
        b = 5,
        l = 5
      )
    )


  ggsave(
    file.path(
      out_dir,
      "SUPP_GSEA_KEGG_peroxisome_curve.png"
    ),
    plot = p_perox,
    width = 5.2,
    height = 4.2,
    dpi = 600,
    bg = "white"
  )


  ggsave(
    file.path(
      out_dir,
      "SUPP_GSEA_KEGG_peroxisome_curve.pdf"
    ),
    plot = p_perox,
    width = 5.2,
    height = 4.2,
    bg = "white"
  )

} else {

  message(
    "Peroxisome was not found in KEGG GSEA results."
  )
}


# ============================================================
# 36. SUMMARY
# ============================================================

summary_gsea <- tibble(

  metric = c(

    "Ranked genes",

    "GO:BP pathways tested",
    "GO:BP significant FDR < 0.25",
    "GO:BP significant UP",
    "GO:BP significant DOWN",

    "GO:CC pathways tested",
    "GO:CC significant FDR < 0.25",
    "GO:CC significant UP",
    "GO:CC significant DOWN",

    "GO:MF pathways tested",
    "GO:MF significant FDR < 0.25",
    "GO:MF significant UP",
    "GO:MF significant DOWN",

    "Mitochondrial GO:CC terms",
    "Significant mitochondrial GO:CC terms",

    "Significant metabolic GO:BP pathways",

    "KEGG ORA DEGs",
    "KEGG ORA background genes",
    "KEGG ORA pathways tested",
    "KEGG ORA significant FDR < 0.05",
    "KEGG ORA Peroxisome detected",

    "Supplementary KEGG GSEA pathways tested",
    "Supplementary KEGG GSEA significant FDR < 0.25",
    "Supplementary KEGG GSEA UP",
    "Supplementary KEGG GSEA DOWN",
    "Supplementary KEGG GSEA Peroxisome detected"
  ),

  value = c(

    length(gene_list),

    nrow(gsea_bp_df),
    nrow(gsea_bp_sig),
    nrow(bp_up),
    nrow(bp_down),

    nrow(gsea_cc_df),
    nrow(gsea_cc_sig),
    nrow(cc_up),
    nrow(cc_down),

    nrow(gsea_mf_df),
    nrow(gsea_mf_sig),
    nrow(mf_up),
    nrow(mf_down),

    nrow(mito_cc),
    nrow(mito_cc_sig),

    nrow(gsea_metabolic_sig),

    nrow(kegg_deg),
    length(kegg_background),
    nrow(kegg_ora_df),
    nrow(kegg_ora_sig),
    nrow(kegg_ora_peroxisome),

    nrow(gsea_kegg_df),
    nrow(gsea_kegg_sig),
    nrow(kegg_up),
    nrow(kegg_down),
    nrow(peroxisome_hits)
  )
)

write_csv(
  summary_gsea,
  file.path(
    out_dir,
    "00_GSEA_KEGG_summary.csv"
  )
)

cat(
  "\n==============================\n"
)

cat(
  "GSEA / KEGG SUMMARY\n"
)

cat(
  "==============================\n"
)

print(
  summary_gsea
)

cat(
  "\nAnalysis completed successfully.\n"
)
