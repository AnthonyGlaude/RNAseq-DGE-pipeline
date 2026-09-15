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
  library(scales)
  library(tidytext)
  library(svglite)
})

# 0. PATHS

stats_file <- snakemake@input[["stats"]]
out_dir <- snakemake@output[["outdir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. PARAMETERS

deg_fdr_cutoff <- 0.05
ora_fdr_cutoff <- 0.05

gsea_fdr_cutoff <- 0.25
min_gs_size <- 15
max_gs_size <- 500

n_top_bp_up   <- 10
n_top_bp_down <- 10
n_top_cc_up   <- 10
n_top_cc_down <- 10
n_top_kegg    <- 10

fdr_colours <- c( "#081D58", "#253494", "#6A1B9A", "#B2182B" )

# 2. READ DESEQ2 RESULTS

if (!file.exists(stats_file)) {
  stop( paste0( "DESeq2 file not found:\n", stats_file ) )
}

df <- read_csv( stats_file, show_col_types = FALSE )

cat("\nColumns in DESeq2 file:\n")
print(names(df))

# 3. DETECT ENSEMBL COLUMN

candidate_gene_columns <- c( "gene", "Gene", "ENSEMBL", "ensembl", "ensembl_gene_id", "gene_id", "...1" )

gene_col <- candidate_gene_columns[ candidate_gene_columns %in% names(df) ][1]

if (is.na(gene_col)) {

  ensembl_match <- vapply(
    df,
    function(x) {
      x <- as.character(x)
      any( grepl("^ENSMUSG", x), na.rm = TRUE )
    },
    logical(1)
  )

  matching_columns <- names(df)[ensembl_match]

  if (length(matching_columns) == 0) {
    stop("Could not identify an Ensembl gene column.")
  }

  gene_col <- matching_columns[1]
}

cat( "\nEnsembl column detected: ", gene_col, "\n", sep = "" )

# 4. CHECK REQUIRED COLUMNS

required_columns <- c( "log2FoldChange", "stat", "padj" )

missing_columns <- setdiff( required_columns, names(df) )

if (length(missing_columns) > 0) {
  stop( paste0( "Missing DESeq2 columns: ", paste(missing_columns, collapse = ", ") ) )
}

# 5. CLEAN IDS + MAP SYMBOL / ENTREZ

df <- df %>%
  mutate( ENSEMBL = sub( "\\..*$", "", as.character(.data[[gene_col]]) ) )

df$SYMBOL <- suppressMessages( AnnotationDbi::mapIds( org.Mm.eg.db, keys = df$ENSEMBL, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first" ) )

df$ENTREZID <- suppressMessages(
  AnnotationDbi::mapIds( org.Mm.eg.db, keys = df$ENSEMBL, keytype = "ENSEMBL", column = "ENTREZID", multiVals = "first" )
)

df$SYMBOL   <- as.character(df$SYMBOL)
df$ENTREZID <- as.character(df$ENTREZID)

cat( "\nGenes with SYMBOL: ", sum(!is.na(df$SYMBOL)), " / ", nrow(df), "\n", sep = "" )

cat( "Genes with ENTREZ: ", sum(!is.na(df$ENTREZID)), " / ", nrow(df), "\n", sep = "" )

# 6. DEFINE DEGs

deg_all <- df %>%
  filter( !is.na(padj), padj < deg_fdr_cutoff, !is.na(ENTREZID), ENTREZID != "" ) %>%
  distinct( ENTREZID, .keep_all = TRUE )

deg_up <- deg_all %>%
  filter(log2FoldChange > 0)

deg_down <- deg_all %>%
  filter(log2FoldChange < 0)

cat( "\nDEGs padj < ", deg_fdr_cutoff, ": ", nrow(deg_all), "\n", sep = "" )

cat( "UP: ", nrow(deg_up), " | DOWN: ", nrow(deg_down), "\n", sep = "" )

write_csv( deg_all, file.path(out_dir, "01_DEGs_FDR005_all.csv") )

write_csv( deg_up, file.path(out_dir, "02_DEGs_FDR005_UP.csv") )

write_csv( deg_down, file.path(out_dir, "03_DEGs_FDR005_DOWN.csv") )

# 7. ORA BACKGROUND

ora_background <- df %>%
  filter( !is.na(padj), !is.na(ENTREZID), ENTREZID != "" ) %>%
  distinct(ENTREZID) %>%
  pull(ENTREZID)

cat( "ORA background: ", length(ora_background), "\n", sep = "" )

# 8. GO ORA FUNCTION

run_go_ora <- function(entrez_ids, ontology) {

  if (length(entrez_ids) == 0) {
    return(NULL)
  }

  enrichGO(
    gene = entrez_ids,
    universe = ora_background,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    readable = TRUE
  )
}

safe_enrich_df <- function(x) {
  if (is.null(x)) {
    return(tibble())
  }
  out <- as.data.frame(x)
  if (nrow(out) == 0) {
    return(tibble())
  }
  as_tibble(out)
}

# 9. RUN GO:BP ORA

cat("\nRunning GO:BP ORA UP...\n")
go_bp_up <- run_go_ora(deg_up$ENTREZID, "BP")

cat("Running GO:BP ORA DOWN...\n")
go_bp_down <- run_go_ora(deg_down$ENTREZID, "BP")

go_bp_up_df <- safe_enrich_df(go_bp_up) %>%
  arrange(p.adjust)

go_bp_down_df <- safe_enrich_df(go_bp_down) %>%
  arrange(p.adjust)

write_csv( go_bp_up_df, file.path(out_dir, "04_GO_BP_ORA_UP_all.csv") )

write_csv( go_bp_down_df, file.path(out_dir, "05_GO_BP_ORA_DOWN_all.csv") )

go_bp_up_sig <- go_bp_up_df %>%
  filter( !is.na(p.adjust), p.adjust < ora_fdr_cutoff )

go_bp_down_sig <- go_bp_down_df %>%
  filter( !is.na(p.adjust), p.adjust < ora_fdr_cutoff )

write_csv( go_bp_up_sig, file.path(out_dir, "06_GO_BP_ORA_UP_FDR005.csv") )

write_csv( go_bp_down_sig, file.path(out_dir, "07_GO_BP_ORA_DOWN_FDR005.csv") )

# 10. GO:BP ORA DOTPLOT FUNCTION -- ONE DIRECTION PER FIGURE

make_ora_single_dotplot <- function(
  enrich_df,
  n_top,
  title,
  filename,
  highlight_lipid = FALSE,
  highlight_amino = FALSE,
  highlight_mito = FALSE
) {

  top_df <- enrich_df %>%
    filter( !is.na(p.adjust) ) %>%
    arrange( p.adjust ) %>%
    slice_head( n = n_top )

  if (nrow(top_df) == 0) {
    message("No pathways available for: ", title)
    return(NULL)
  }

  top_df <- top_df %>%
    mutate(

      GeneRatio_numeric =
        as.numeric( sub( "/.*", "", GeneRatio ) ) /
        as.numeric( sub( ".*/", "", GeneRatio ) ),

      is_lipid = str_detect(
        Description,
        regex( paste0( "fatty acid|", "lipid|", "lipoprotein|", "acyl|", "cholesterol|", "triglycer|", "phospholipid|", "PPAR" ), ignore_case = TRUE )
      ),

      is_amino = str_detect(
        Description,
        regex( paste0( "amino acid|", "valine|", "leucine|", "isoleucine|", "alanine|", "aspartate|", "glutamate|", "branched" ), ignore_case = TRUE )
      ),

      is_mito = str_detect( Description, regex( "mitochond", ignore_case = TRUE ) ),

      Description_short = str_wrap( Description, width = 42 ),

      Description_label = case_when(

        highlight_lipid & is_lipid ~ paste0( "<span style='color:#D55E00;font-weight:bold'>", Description_short, "</span>" ),

        highlight_amino & is_amino ~ paste0( "<span style='color:#0072B2;font-weight:bold'>", Description_short, "</span>" ),

        highlight_mito & is_mito ~ paste0( "<span style='color:#7B3294;font-weight:bold'>", Description_short, "</span>" ),

        TRUE ~ Description_short
      )
    )

  top_df <- top_df %>%
    arrange( GeneRatio_numeric ) %>%
    mutate( Description_plot = factor( Description_label, levels = Description_label ) )

  p <- ggplot( top_df, aes( x = GeneRatio_numeric, y = Description_plot ) ) +

    geom_point( aes( size = Count, colour = p.adjust ) ) +

    scale_colour_gradientn( colours = fdr_colours, trans = "log10" ) +

    scale_size_continuous( range = c(2.4, 5.8), breaks = pretty( top_df$Count, n = 3 ) ) +

    scale_x_continuous(
      breaks = pretty( range( top_df$GeneRatio_numeric, na.rm = TRUE ), n = 5 ),
      labels = scales::label_number( accuracy = 0.01 ),
      expand = expansion( mult = c(0.03, 0.08) )
    ) +

    theme_classic( base_size = 10 ) +

    theme(
      axis.text.y = ggtext::element_markdown( size = 11, lineheight = 0.98 ),
      axis.text.x = element_text( size = 8.5 ),
      axis.title.x = element_text( size = 10 ),
      plot.title = element_text( hjust = 0.5, size = 11.5, face = "bold" ),
      legend.title = element_text( size = 9 ),
      legend.text = element_text( size = 8 ),
      legend.position = "right",
      plot.margin = margin( t = 8, r = 8, b = 6, l = 6 )
    ) +

    labs( title = title, x = "GeneRatio", y = NULL, size = "Count", colour = "FDR" )

  ggsave( file.path( out_dir, paste0(filename, ".png") ), p, width = 7.4, height = 4.8, dpi = 600, bg = "white" )

  ggsave( file.path( out_dir, paste0(filename, ".pdf") ), p, width = 7.4, height = 4.8, bg = "white" )

  ggsave( file.path( out_dir, paste0(filename, ".svg") ), p, width = 7.4, height = 4.8, bg = "white" )

  write_csv(
    top_df %>%
      select( ID, Description, GeneRatio, GeneRatio_numeric, BgRatio, Count, pvalue, p.adjust ),
    file.path( out_dir, paste0(filename, "_table.csv") )
  )

  return(p)
}

# 11. PANELS C AND D -- GO:BP ORA

p_bp_up <- make_ora_single_dotplot(
  enrich_df = go_bp_up_sig,
  n_top = n_top_bp_up,
  title = "Top pathways enriched in upregulated genes in HFpEF vs. Chow",
  filename = "08_PANEL_C_GO_BP_ORA_UP_top10",
  highlight_lipid = TRUE,
  highlight_amino = TRUE,
  highlight_mito = FALSE
)

p_bp_down <- make_ora_single_dotplot(
  enrich_df = go_bp_down_sig,
  n_top = n_top_bp_down,
  title = "Top pathways enriched in downregulated genes in HFpEF vs. Chow",
  filename = "09_PANEL_D_GO_BP_ORA_DOWN_top10",
  highlight_lipid = TRUE,
  highlight_amino = TRUE,
  highlight_mito = TRUE
)

# 11B. COMBINED UP/DOWN DOTPLOT FUNCTION FOR GO:CC

make_ora_updown_dotplot <- function(
  up_df,
  down_df,
  n_up,
  n_down,
  title,
  filename,
  highlight_lipid = FALSE,
  highlight_amino = FALSE,
  highlight_mito = FALSE,
  panel_labels = c("Up-regulated", "Down-regulated")
) {

  # 1. TOP pathways separately by FDR

  top_up <- up_df %>%
    arrange(p.adjust) %>%
    slice_head(n = n_up) %>%
    mutate( direction = panel_labels[1] )

  top_down <- down_df %>%
    arrange(p.adjust) %>%
    slice_head(n = n_down) %>%
    mutate( direction = panel_labels[2] )

  # 2. Combine

  top_df <- bind_rows( top_up, top_down )

  if (nrow(top_df) == 0) {

    message( "No pathways available for: ", title )

    return(NULL)
  }

  # 3. GeneRatio + pathway highlighting

  top_df <- top_df %>%

    mutate(

      GeneRatio_numeric =
        as.numeric( sub( "/.*", "", GeneRatio ) ) /
        as.numeric( sub( ".*/", "", GeneRatio ) ),

      is_lipid = str_detect( Description, regex( "fatty acid|lipid|lipoprotein|acyl|cholesterol|triglycer|phospholipid|PPAR", ignore_case = TRUE ) ),

      is_amino = str_detect( Description, regex( "amino acid|valine|leucine|isoleucine|alanine|aspartate|glutamate|branched", ignore_case = TRUE ) ),

      is_mito = str_detect( Description, regex( "mitochond", ignore_case = TRUE ) ),

      Description_label = case_when(

        highlight_lipid & is_lipid ~ paste0( "<span style='color:#D55E00'>", Description, "</span>" ),

        highlight_amino & is_amino ~ paste0( "<span style='color:#0072B2'>", Description, "</span>" ),

        highlight_mito & is_mito ~ paste0( "<span style='color:#0072B2'>", Description, "</span>" ),

        TRUE ~ Description
      )
    )

  # 4. Order pathways globally

  pathway_order <- top_df %>%

    group_by( Description_label ) %>%

    summarise( best_fdr = min( p.adjust, na.rm = TRUE ), .groups = "drop" ) %>%

    arrange( desc(best_fdr) ) %>%

    pull( Description_label )

  top_df <- top_df %>%

    mutate(

      Description_plot = factor( Description_label, levels = pathway_order ),

      direction = factor( direction, levels = panel_labels )
    )

  # 5. Plot

  p <- ggplot(

    top_df,

    aes( x = GeneRatio_numeric, y = Description_plot )
  ) +

    geom_point(

      aes( size = Count, colour = p.adjust )
    ) +

    facet_grid( . ~ direction ) +

    scale_colour_gradientn(

      colours = c( "#081D58", "#253494", "#6A1B9A", "#B2182B" ),

      trans = "log10"
    ) +

    scale_size_continuous(

      breaks = pretty( top_df$Count, n = 4 ),

      range = c( 2.5, 8 )
    ) +

    scale_x_continuous(

      breaks = pretty( range( top_df$GeneRatio_numeric, na.rm = TRUE ), n = 5 ),

      labels = scales::label_number( accuracy = 0.01 ),

      expand = expansion( mult = c( 0.08, 0.12 ) )
    ) +

    theme_bw( base_size = 11 ) +

    theme(

      panel.grid.major.y = element_blank(),

      panel.grid.minor = element_blank(),

      panel.grid.major.x = element_line( colour = "grey90", linewidth = 0.4 ),

      axis.text.y = ggtext::element_markdown( size = 8.5 ),

      axis.text.x = element_text( size = 9 ),

      axis.title.x = element_text( size = 11 ),

      strip.background = element_rect( fill = "white", colour = "grey60" ),

      strip.text = element_text( size = 10, face = "bold" ),

      plot.title = element_text( hjust = 0.5, size = 13 ),

      legend.position = "right",

      plot.margin = margin( t = 10, r = 15, b = 10, l = 10 )
    ) +

    labs(

      title = title,

      x = "GeneRatio",

      y = NULL,

      size = "Count",

      colour = "FDR"
    )

  # 6. SAVE

  ggsave(

    file.path( out_dir, paste0( filename, ".png" ) ),

    p,

    width = 10,

    height = 7,

    dpi = 600,

    bg = "white"
  )

  ggsave(

    file.path( out_dir, paste0( filename, ".pdf" ) ),

    p,

    width = 10,

    height = 7,

    bg = "white"
  )

  ggsave(

    file.path( out_dir, paste0( filename, ".svg" ) ),

    p,

    width = 10,

    height = 7,

    bg = "white"
  )

  return(p)
}

# 12. RUN GO:CC ORA

cat("\nRunning GO:CC ORA UP...\n")
go_cc_up <- run_go_ora(deg_up$ENTREZID, "CC")

cat("Running GO:CC ORA DOWN...\n")
go_cc_down <- run_go_ora(deg_down$ENTREZID, "CC")

go_cc_up_df <- safe_enrich_df(go_cc_up) %>%
  arrange(p.adjust)

go_cc_down_df <- safe_enrich_df(go_cc_down) %>%
  arrange(p.adjust)

write_csv( go_cc_up_df, file.path(out_dir, "10_GO_CC_ORA_UP_all.csv") )

write_csv( go_cc_down_df, file.path(out_dir, "11_GO_CC_ORA_DOWN_all.csv") )

go_cc_up_sig <- go_cc_up_df %>%
  filter( !is.na(p.adjust), p.adjust < ora_fdr_cutoff )

go_cc_down_sig <- go_cc_down_df %>%
  filter( !is.na(p.adjust), p.adjust < ora_fdr_cutoff )

write_csv( go_cc_up_sig, file.path(out_dir, "12_GO_CC_ORA_UP_FDR005.csv") )

write_csv( go_cc_down_sig, file.path(out_dir, "13_GO_CC_ORA_DOWN_FDR005.csv") )

# 13. PANEL E -- GO:CC ORA DOTPLOT

p_cc_ora <- make_ora_updown_dotplot(
  up_df = go_cc_up_sig,
  down_df = go_cc_down_sig,
  n_up = n_top_cc_up,
  n_down = n_top_cc_down,
  title = "GO enrichment analysis (CC)",
  filename = "14_PANEL_F_GO_CC_ORA_top10_UP_DOWN",
  highlight_lipid = FALSE,
  highlight_amino = FALSE,
  highlight_mito = TRUE,
  panel_labels = c("Up-regulated", "Down-regulated")
)

# 14. TARGETED GSEA FOR A BIOLOGICALLY SELECTED GO:BP PATHWAY

rank_df <- df %>%
  filter( !is.na(SYMBOL), SYMBOL != "", !is.na(stat) ) %>%
  mutate(rank_score = stat) %>%
  arrange(desc(abs(rank_score))) %>%
  distinct(SYMBOL, .keep_all = TRUE)

gene_list <- rank_df$rank_score
names(gene_list) <- rank_df$SYMBOL
gene_list <- sort(gene_list, decreasing = TRUE)

n_ties <- sum(duplicated(gene_list))
if (n_ties > 0) {
  gene_list <- gene_list + seq_along(gene_list) * 1e-12
  gene_list <- sort(gene_list, decreasing = TRUE)
}

write_csv( tibble( SYMBOL = names(gene_list), rank_score = as.numeric(gene_list) ), file.path(out_dir, "15_GSEA_targeted_preranked_gene_list.csv") )

set.seed(12345)

gsea_bp <- clusterProfiler::gseGO(
  geneList = gene_list,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = min_gs_size,
  maxGSSize = max_gs_size,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  seed = TRUE,
  by = "fgsea"
)

gsea_bp_df <- as.data.frame(gsea_bp)

write_csv( gsea_bp_df, file.path(out_dir, "16_GSEA_BP_targeted_all.csv") )

fatty_hits <- gsea_bp_df %>%
  filter( str_detect( Description, regex( "^fatty acid metabolic process$", ignore_case = TRUE ) ) ) %>%
  arrange(p.adjust)

write_csv( fatty_hits, file.path(out_dir, "17_GSEA_BP_fatty_acid_metabolic_process.csv") )

cat( "\nFatty acid metabolic process hits: ", nrow(fatty_hits), "\n", sep = "" )

if (nrow(fatty_hits) > 0) {

  print( fatty_hits %>% dplyr::select( ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust ) )

  fatty_id <- fatty_hits$ID[1]
  fatty_nes <- round(fatty_hits$NES[1], 2)
  fatty_fdr <- formatC( fatty_hits$p.adjust[1], format = "g", digits = 3 )

  p_fatty <- enrichplot::gseaplot2(
    gsea_bp,
    geneSetID = fatty_id,
    title = "[GO:BP] Fatty acid metabolic process",
    pvalue_table = FALSE,
    ES_geom = "line",
    subplots = c(1, 2),
    rel_heights = c(1.8, 0.45),
    base_size = 11
  )

  clean_panel <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect( colour = "black", fill = NA, linewidth = 0.6 ),
    axis.line = element_blank(),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.background = element_rect(fill = "white", colour = NA)
  )

  p_fatty[[1]] <- p_fatty[[1]] +
    labs(y = "Enrichment Score") +
    clean_panel +
    theme(
      plot.title = element_text( hjust = 0.5, size = 12, face = "plain" ),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
    ) +
    annotate( "text", x = Inf, y = Inf, label = paste0( "NES = ", fatty_nes, "\nFDR = ", fatty_fdr ), hjust = 1.08, vjust = 1.15, size = 3.8 )

  p_fatty[[2]] <- p_fatty[[2]] +
    labs( x = "Rank in Ordered Dataset" ) +
    clean_panel +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
    )

  ggsave( file.path(out_dir, "18_PANEL_E_GSEA_fatty_acid_metabolic_process.png"), plot = p_fatty, width = 5.2, height = 4.2, dpi = 600, bg = "white" )

  ggsave( file.path(out_dir, "18_PANEL_E_GSEA_fatty_acid_metabolic_process.pdf"), plot = p_fatty, width = 5.2, height = 4.2, bg = "white" )

  ggsave( file.path(out_dir, "18_PANEL_E_GSEA_fatty_acid_metabolic_process.svg"), plot = p_fatty, width = 5.2, height = 3.5, bg = "white" )

} else {
  warning("Fatty acid metabolic process was not found in gseGO results.")
}

# 15. KEGG ORA GLOBAL

cat("\nRunning KEGG ORA...\n")

kegg_ora <- clusterProfiler::enrichKEGG(
  gene = deg_all$ENTREZID,
  universe = ora_background,
  organism = "mmu",
  keyType = "ncbi-geneid",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1
)

kegg_ora_df <- safe_enrich_df(kegg_ora) %>%
  arrange(p.adjust)

write_csv( kegg_ora_df, file.path(out_dir, "19_KEGG_ORA_all.csv") )

kegg_ora_sig <- kegg_ora_df %>%
  filter( !is.na(p.adjust), p.adjust < ora_fdr_cutoff )

write_csv( kegg_ora_sig, file.path(out_dir, "20_KEGG_ORA_FDR005.csv") )

# 16. PANEL F -- KEGG ORA DOTPLOT

if (nrow(kegg_ora_df) > 0) {

  kegg_top <- kegg_ora_df %>%
    filter( !is.na(p.adjust) ) %>%
    arrange( p.adjust ) %>%
    slice_head( n = n_top_kegg ) %>%
    mutate(

      GeneRatio_numeric =
        as.numeric( sub( "/.*", "", GeneRatio ) ) /
        as.numeric( sub( ".*/", "", GeneRatio ) ),

      is_lipid = str_detect(
        Description,
        regex( paste0( "fatty acid|", "lipid|", "PPAR|", "peroxisome|", "acyl|", "cholesterol|", "triglycer|", "phospholipid" ), ignore_case = TRUE )
      ),

      is_amino = str_detect(
        Description,
        regex(
          paste0( "amino acid|", "valine|", "leucine|", "isoleucine|", "tryptophan|", "alanine|", "glutamate|", "aspartate" ),
          ignore_case = TRUE
        )
      ),

      is_apelin = str_detect( Description, regex( "apelin", ignore_case = TRUE ) ),

      Description_short = str_wrap( Description, width = 38 ),

      Description_label = case_when(

        is_apelin ~ paste0( "<span style='color:#009E73;font-weight:bold'>", Description_short, "</span>" ),

        is_lipid ~ paste0( "<span style='color:#D55E00;font-weight:bold'>", Description_short, "</span>" ),

        is_amino ~ paste0( "<span style='color:#0072B2;font-weight:bold'>", Description_short, "</span>" ),

        TRUE ~ Description_short
      )
    )

  # IMPORTANT:

  kegg_top <- kegg_top %>%
    arrange( GeneRatio_numeric ) %>%
    mutate( Description_plot = factor( Description_label, levels = Description_label ) )

  p_kegg_ora <- ggplot( kegg_top, aes( x = GeneRatio_numeric, y = Description_plot ) ) +

    geom_point( aes( size = Count, colour = p.adjust ) ) +

    scale_colour_gradientn( colours = fdr_colours, trans = "log10" ) +

    scale_size_continuous( range = c( 2.2, 5.5 ), breaks = pretty( kegg_top$Count, n = 3 ) ) +

    scale_x_continuous(
      breaks = pretty( range( kegg_top$GeneRatio_numeric, na.rm = TRUE ), n = 5 ),
      labels = scales::label_number( accuracy = 0.01 ),
      expand = expansion( mult = c( 0.03, 0.08 ) )
    ) +

    theme_classic( base_size = 10 ) +

    theme(

      axis.text.y = ggtext::element_markdown( size = 10.8, lineheight = 0.98 ),

      axis.text.x = element_text( size = 8 ),

      axis.title.x = element_text( size = 10 ),

      plot.title = element_text( hjust = 0.5, size = 12, face = "bold" ),

      legend.title = element_text( size = 9 ),

      legend.text = element_text( size = 8 ),

      legend.position = "right",

      plot.margin = margin( t = 8, r = 8, b = 6, l = 6 )
    ) +

    labs( title = "KEGG pathway enrichment", x = "GeneRatio", y = NULL, size = "Gene count", colour = "FDR" )

  ggsave( file.path( out_dir, "21_PANEL_G_KEGG_ORA_top10_dotplot.png" ), p_kegg_ora, width = 8.2, height = 4.8, dpi = 600, bg = "white" )

  ggsave( file.path( out_dir, "21_PANEL_G_KEGG_ORA_top10_dotplot.pdf" ), p_kegg_ora, width = 8.2, height = 4.8, bg = "white" )

  ggsave( file.path( out_dir, "21_PANEL_G_KEGG_ORA_top10_dotplot.svg" ), p_kegg_ora, width = 8.2, height = 4.8, bg = "white" )

  write_csv(
    kegg_top %>%
      select( ID, Description, GeneRatio, GeneRatio_numeric, BgRatio, Count, pvalue, p.adjust ),
    file.path( out_dir, "21_PANEL_G_KEGG_ORA_top10_table.csv" )
  )
}

# 17. REPORT SPECIFIC KEGG TERMS

kegg_peroxisome <- kegg_ora_df %>%
  filter( ID == "mmu04146" | str_detect( Description, regex( "^Peroxisome$", ignore_case = TRUE ) ) )

kegg_ppar <- kegg_ora_df %>%
  filter( str_detect( Description, regex( "^PPAR signaling pathway$", ignore_case = TRUE ) ) )

kegg_fatty <- kegg_ora_df %>%
  filter( str_detect( Description, regex( "fatty acid", ignore_case = TRUE ) ) )

kegg_apelin <- kegg_ora_df %>%
  filter( str_detect( Description, regex( "^Apelin signaling pathway$", ignore_case = TRUE ) ) )

write_csv( kegg_peroxisome, file.path( out_dir, "22_KEGG_ORA_peroxisome.csv" ) )

write_csv( kegg_ppar, file.path( out_dir, "23_KEGG_ORA_PPAR_signaling.csv" ) )

write_csv( kegg_fatty, file.path( out_dir, "24_KEGG_ORA_fatty_acid_pathways.csv" ) )

write_csv( kegg_apelin, file.path( out_dir, "25_KEGG_ORA_apelin_signaling_pathway.csv" ) )

if (nrow(kegg_peroxisome) > 0) {

  cat( "\nKEGG ORA -- PEROXISOME\n" )

  print( as_tibble( kegg_peroxisome %>% dplyr::select( ID, Description, GeneRatio, BgRatio, Count, pvalue, p.adjust ) ) )
}

if (nrow(kegg_ppar) > 0) {

  cat( "\nKEGG ORA -- PPAR SIGNALING\n" )

  print( as_tibble( kegg_ppar %>% dplyr::select( ID, Description, GeneRatio, BgRatio, Count, pvalue, p.adjust ) ) )
}

if (nrow(kegg_apelin) > 0) {

  cat( "\nKEGG ORA -- APELIN SIGNALING PATHWAY\n" )

  print( as_tibble( kegg_apelin %>% dplyr::select( ID, Description, GeneRatio, BgRatio, Count, pvalue, p.adjust ) ) )
}

# 17B. GSEA KEGG -- PEROXISOME / PPAR

rank_df_kegg <- df %>%
  filter( !is.na(ENTREZID), ENTREZID != "", !is.na(stat) ) %>%
  mutate(rank_score = stat) %>%
  arrange(desc(abs(rank_score))) %>%
  distinct(ENTREZID, .keep_all = TRUE)

gene_list_kegg <- rank_df_kegg$rank_score
names(gene_list_kegg) <- rank_df_kegg$ENTREZID
gene_list_kegg <- sort(gene_list_kegg, decreasing = TRUE)

n_ties_kegg <- sum(duplicated(gene_list_kegg))

if (n_ties_kegg > 0) {
  gene_list_kegg <- gene_list_kegg +
    seq_along(gene_list_kegg) * 1e-12
  gene_list_kegg <- sort( gene_list_kegg, decreasing = TRUE )
}

write_csv(
  tibble( ENTREZID = names(gene_list_kegg), rank_score = as.numeric(gene_list_kegg) ),
  file.path( out_dir, "26_GSEA_KEGG_preranked_gene_list.csv" )
)

cat("\nRunning KEGG GSEA for targeted pathway plots...\n")

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

gsea_kegg_df <- safe_enrich_df( gsea_kegg ) %>%
  arrange(p.adjust)

write_csv( gsea_kegg_df, file.path( out_dir, "27_GSEA_KEGG_all.csv" ) )

make_targeted_kegg_gsea <- function(
  pathway_regex,
  plot_title,
  filename
) {

  pathway_hit <- gsea_kegg_df %>%
    filter( str_detect( Description, regex( pathway_regex, ignore_case = TRUE ) ) ) %>%
    arrange(p.adjust)

  write_csv( pathway_hit, file.path( out_dir, paste0(filename, "_table.csv") ) )

  if (nrow(pathway_hit) == 0) {
    warning( plot_title, " was not found in KEGG GSEA results." )
    return(NULL)
  }

  pathway_id <- pathway_hit$ID[1]

  pathway_nes <- round( pathway_hit$NES[1], 2 )

  pathway_fdr <- formatC( pathway_hit$p.adjust[1], format = "g", digits = 3 )

  p_pathway <- enrichplot::gseaplot2(
    gsea_kegg,
    geneSetID = pathway_id,
    title = plot_title,
    pvalue_table = FALSE,
    ES_geom = "line",
    subplots = c(1, 2),
    rel_heights = c(1.8, 0.45),
    base_size = 11
  )

  pathway_clean_panel <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect( colour = "black", fill = NA, linewidth = 0.6 ),
    axis.line = element_blank(),
    axis.text = element_text( size = 9 ),
    axis.title = element_text( size = 10 ),
    plot.background = element_rect( fill = "white", colour = NA )
  )

  p_pathway[[1]] <- p_pathway[[1]] +
    labs( y = "Enrichment Score" ) +
    pathway_clean_panel +
    theme(
      plot.title = element_text( hjust = 0.5, size = 12, face = "plain" ),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin( t = 5, r = 5, b = 0, l = 5 )
    ) +
    annotate( "text", x = Inf, y = Inf, label = paste0( "NES = ", pathway_nes, "\nFDR = ", pathway_fdr ), hjust = 1.08, vjust = 1.15, size = 3.8 )

  p_pathway[[2]] <- p_pathway[[2]] +
    labs( x = "Rank in Ordered Dataset" ) +
    pathway_clean_panel +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin( t = 0, r = 5, b = 5, l = 5 )
    )

  ggsave( file.path( out_dir, paste0(filename, ".png") ), plot = p_pathway, width = 5.4, height = 3.5, dpi = 600, bg = "white" )

  ggsave( file.path( out_dir, paste0(filename, ".pdf") ), plot = p_pathway, width = 5.4, height = 3.5, bg = "white" )

  ggsave( file.path( out_dir, paste0(filename, ".svg") ), plot = p_pathway, width = 5.4, height = 3.5, bg = "white" )

  return(p_pathway)
}

p_gsea_peroxisome <- make_targeted_kegg_gsea( pathway_regex = "^Peroxisome$", plot_title = "[KEGG] Peroxisome", filename = "28_GSEA_KEGG_peroxisome" )

p_gsea_ppar <- make_targeted_kegg_gsea(
  pathway_regex = "^PPAR signaling pathway$",
  plot_title = "[KEGG] PPAR signaling pathway",
  filename = "29_GSEA_KEGG_PPAR_signaling_pathway"
)

# 17C. PATHVIEW -- APELIN SIGNALING PATHWAY
# Purpose:

pathview_dir <- file.path( out_dir, '30_pathview_apelin_signaling' )

dir.create( pathview_dir, recursive = TRUE, showWarnings = FALSE )

if (nrow(kegg_apelin) > 0 && 'geneID' %in% colnames(kegg_apelin)) {

  apelin_entrez <- unique( unlist( strsplit( kegg_apelin$geneID[1], '/', fixed = TRUE ) ) )

  apelin_gene_table <- df %>%
    filter( !is.na(ENTREZID), ENTREZID %in% apelin_entrez ) %>%
    select( ENTREZID, SYMBOL, ENSEMBL, log2FoldChange, pvalue, padj ) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    distinct(ENTREZID, .keep_all = TRUE)

  write_csv( apelin_gene_table, file.path( pathview_dir, '30A_KEGG_Apelin_signaling_gene_members.csv' ) )
}

pathview_df <- df %>%
  filter( !is.na(ENTREZID), ENTREZID != '', !is.na(log2FoldChange) ) %>%
  mutate( abs_log2FC = abs(log2FoldChange) ) %>%
  arrange(desc(abs_log2FC)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

pathview_gene_data <- pathview_df$log2FoldChange
names(pathview_gene_data) <- pathview_df$ENTREZID

pathview_limit <- max( 1, unname( quantile( abs(pathview_gene_data), probs = 0.95, na.rm = TRUE ) ) )

if (nrow(kegg_apelin) > 0) {

  if (requireNamespace('pathview', quietly = TRUE)) {

    cat('\nRunning Pathview for KEGG Apelin signaling pathway...\n')

    old_wd <- getwd()
    setwd(pathview_dir)

    pv_out <- tryCatch(
      {
        pathview::pathview(
          gene.data  = pathview_gene_data,
          pathway.id = '04371',
          species    = 'mmu',
          gene.idtype = 'entrez',
          kegg.native = TRUE,
          same.layer = TRUE,
          out.suffix = 'HFpEF_WT_vs_Chow_Apelin',
          limit = list(gene = pathview_limit, cpd = 1),
          low = list(gene = '#2C7BB6', cpd = 'white'),
          mid = list(gene = 'grey90', cpd = 'white'),
          high = list(gene = '#D7191C', cpd = 'white')
        )
      },
      error = function(e) {
        message('Pathview failed: ', conditionMessage(e))
        NULL
      },
      finally = {
        setwd(old_wd)
      }
    )

    saveRDS( pv_out, file.path( pathview_dir, '30B_pathview_apelin_result.rds' ) )

    generated_png <- list.files( pathview_dir, pattern = '^mmu04371\\..*Apelin.*\\.png$', full.names = TRUE )

    generated_xml <- list.files( pathview_dir, pattern = '^mmu04371\\.xml$', full.names = TRUE )

    generated_pdf <- list.files( pathview_dir, pattern = '^mmu04371\\..*Apelin.*\\.pdf$', full.names = TRUE )

    if (length(generated_png) > 0) {
      file.copy( generated_png[1], file.path( pathview_dir, '30C_PATHVIEW_Apelin_signaling_mmu04371.png' ), overwrite = TRUE )
    }

    if (length(generated_pdf) > 0) {
      file.copy( generated_pdf[1], file.path( pathview_dir, '30D_PATHVIEW_Apelin_signaling_mmu04371.pdf' ), overwrite = TRUE )
    }

    if (length(generated_xml) > 0) {
      file.copy( generated_xml[1], file.path( pathview_dir, '30E_PATHVIEW_Apelin_signaling_mmu04371.xml' ), overwrite = TRUE )
    }

  } else {

    message( '\npathview package not installed. ', 'To enable the Apelin pathway overlay, run:\n', 'BiocManager::install("pathview")\n' )
  }
}

# 18. SUMMARY

summary_enrichment <- tibble(
  metric = c(
    "DEGs padj < 0.05",
    "DEGs UP",
    "DEGs DOWN",
    "ORA background",
    "GO:BP significant UP FDR < 0.05",
    "GO:BP significant DOWN FDR < 0.05",
    "GO:CC significant UP FDR < 0.05",
    "GO:CC significant DOWN FDR < 0.05",
    "Targeted GSEA fatty acid metabolic process hits",
    "KEGG pathways tested",
    "KEGG significant FDR < 0.05",
    "KEGG Peroxisome detected",
    "KEGG PPAR signaling detected",
    "KEGG fatty acid pathways detected",
    "KEGG Apelin signaling pathway detected",
    "Pathview Apelin gene members exported",
    "Pathview Apelin PNG generated"
  ),
  value = c(
    nrow(deg_all),
    nrow(deg_up),
    nrow(deg_down),
    length(ora_background),
    nrow(go_bp_up_sig),
    nrow(go_bp_down_sig),
    nrow(go_cc_up_sig),
    nrow(go_cc_down_sig),
    nrow(fatty_hits),
    nrow(kegg_ora_df),
    nrow(kegg_ora_sig),
    nrow(kegg_peroxisome),
    nrow(kegg_ppar),
    nrow(kegg_fatty),
    nrow(kegg_apelin),
    if (exists("apelin_gene_table")) nrow(apelin_gene_table) else 0,
    as.integer(file.exists(file.path(pathview_dir, "30C_PATHVIEW_Apelin_signaling_mmu04371.png")))
  )
)

write_csv( summary_enrichment, file.path(out_dir, "00_enrichment_summary.csv") )

cat("\n==============================\n")
cat("ENRICHMENT SUMMARY\n")
cat("==============================\n")
print(summary_enrichment)
cat("\nAnalysis completed successfully :D \n")
