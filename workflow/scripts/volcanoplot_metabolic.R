#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

# Inputs and outputs
stats_file <- snakemake@input[["stats"]]
enrichr_up_file <- snakemake@input[["enrichr_up"]]
enrichr_down_file <- snakemake@input[["enrichr_down"]]
volcano_dir <- dirname(snakemake@output[["lipid_png"]])
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters
padj_cutoff <- 0.05
y_cap <- 80
y_cap_position <- 78
y_breaks <- c(0, 20, 40, 60, 80)

# GO terms
go_lipid <- c(
  "GO:0006635", # fatty acid beta-oxidation
  "GO:0009062", # fatty acid catabolic process
  "GO:0015909", # long-chain fatty acid transport
  "GO:0006631", # fatty acid metabolic process
  "GO:0019395", # fatty acid oxidation
  "GO:1902001", # fatty acid transmembrane transport
  "GO:0000038", # very long-chain fatty acid metabolic process
  "GO:0032365"  # intracellular lipid transport
)

go_amino_acid <- c(
  "GO:0009081", # branched-chain amino acid metabolic process
  "GO:0009063", # amino acid catabolic process
  "GO:0009083", # branched-chain amino acid catabolic process
  "GO:0006551", # L-leucine metabolic process
  "GO:1901605", # alpha-amino acid metabolic process
  "GO:0009080", # pyruvate family amino acid catabolic process
  "GO:0006520"  # amino acid metabolic process
)

go_glucose <- c(
  "GO:0006006", # glucose metabolic process
  "GO:0006007", # glucose catabolic process
  "GO:0006096"  # glycolytic process
)

# Read and annotate DESeq2 results
df <- read_csv(stats_file, show_col_types = FALSE)
cat("\nColumns in DESeq2 file:\n")
print(names(df))

gene_candidates <- c("gene", "gene_id", "ensembl_gene_id", "ENSEMBL", "...1")
gene_col <- gene_candidates[gene_candidates %in% names(df)][1]

if (is.na(gene_col)) {
  stop("Could not identify gene ID column. Columns are: ", paste(names(df), collapse = ", "))
}

cat("\nGene ID column detected: ", gene_col, "\n", sep = "")
df <- df %>% mutate(ENSEMBL = sub("\\..*$", "", as.character(.data[[gene_col]])))

df$SYMBOL <- suppressMessages(AnnotationDbi::mapIds(
  org.Mm.eg.db, keys = df$ENSEMBL, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"
))

df <- df %>% mutate(SYMBOL = as.character(SYMBOL), SYMBOL_KEY = toupper(SYMBOL))
cat("\nGenes with SYMBOL: ", sum(!is.na(df$SYMBOL)), " / ", nrow(df), "\n", sep = "")

# Prepare volcano data
res_plot <- df %>%
  filter(!is.na(log2FoldChange), !is.na(padj), padj > 0) %>%
  mutate(
    y_real = -log10(padj),
    above_y_cap = y_real > y_cap,
    y_plot = if_else(above_y_cap, y_cap_position, y_real)
  )

cat("\nGenes included: ", nrow(res_plot), "\n", sep = "")
cat("Genes above Y cap: ", sum(res_plot$above_y_cap), "\n", sep = "")
cat("Maximum -log10(FDR): ", max(res_plot$y_real, na.rm = TRUE), "\n", sep = "")
cat("log2FC range: ", paste(range(res_plot$log2FoldChange, na.rm = TRUE), collapse = " "), "\n", sep = "")

max_abs_x <- max(abs(res_plot$log2FoldChange), na.rm = TRUE)
x_limit <- ceiling(max_abs_x * 1.05)
if (x_limit < 1) x_limit <- 1
cat("Dynamic X limit: ", -x_limit, " to ", x_limit, "\n", sep = "")

# Read Enrichr results
enrichr_up <- read_csv(enrichr_up_file, show_col_types = FALSE)
enrichr_down <- read_csv(enrichr_down_file, show_col_types = FALSE)

find_column <- function(df, candidates) {
  candidate <- candidates[candidates %in% names(df)]
  if (length(candidate) == 0) stop("Could not find column among: ", paste(candidates, collapse = ", "))
  candidate[1]
}

term_col <- find_column(enrichr_up, c("Term", "term", "Description", "description"))
genes_col <- find_column(enrichr_up, c("Genes", "genes", "geneID", "Genes;"))

extract_enrichr_genes <- function(res_table, go_ids, term_col, genes_col) {
  term_vector <- as.character(res_table[[term_col]])
  keep <- vapply(term_vector, function(term) {
    any(vapply(go_ids, function(go) grepl(go, term, fixed = TRUE), logical(1)))
  }, logical(1))

  rows <- res_table[keep, , drop = FALSE]
  if (nrow(rows) == 0) return(character(0))

  genes <- trimws(unlist(strsplit(as.character(rows[[genes_col]]), "[;,/]")))
  unique(genes[!is.na(genes) & genes != ""])
}

# Extract category genes
lipid_up <- extract_enrichr_genes(enrichr_up, go_lipid, term_col, genes_col)
lipid_down <- extract_enrichr_genes(enrichr_down, go_lipid, term_col, genes_col)
lipid_genes <- unique(c(lipid_up, lipid_down))
lipid_keys <- toupper(lipid_genes)

amino_up <- extract_enrichr_genes(enrichr_up, go_amino_acid, term_col, genes_col)
amino_down <- extract_enrichr_genes(enrichr_down, go_amino_acid, term_col, genes_col)
amino_genes <- unique(c(amino_up, amino_down))
amino_keys <- toupper(amino_genes)

glucose_up <- extract_enrichr_genes(enrichr_up, go_glucose, term_col, genes_col)
glucose_down <- extract_enrichr_genes(enrichr_down, go_glucose, term_col, genes_col)
glucose_genes <- unique(c(glucose_up, glucose_down))
glucose_keys <- toupper(glucose_genes)

cat(
  "\n============================================\n",
  "METABOLIC GENE SETS\n",
  "============================================\n",
  "Lipid genes      : ", length(lipid_genes), "\n",
  "Amino acid genes : ", length(amino_genes), "\n",
  "Glucose genes    : ", length(glucose_genes), "\n",
  sep = ""
)

# Category membership; display priority: lipid > amino acid > glucose
plot_metabolism <- res_plot %>%
  mutate(
    lipid = SYMBOL_KEY %in% lipid_keys,
    amino_acid = SYMBOL_KEY %in% amino_keys,
    glucose = SYMBOL_KEY %in% glucose_keys,
    metabolic_class = case_when(
      lipid ~ "Lipid",
      amino_acid ~ "Amino acid",
      glucose ~ "Glucose",
      TRUE ~ "Other"
    )
  )

cat("\nMetabolic classification:\n")
print(table(plot_metabolism$metabolic_class))

apelin_points <- plot_metabolism %>% filter(SYMBOL %in% c("Apln", "Aplnr"))
cat("\nApln / Aplnr:\n")
print(apelin_points %>% dplyr::select(SYMBOL, log2FoldChange, padj, y_real))

make_metabolic_volcano <- function(data, selected_column, output_paths) {
  selected_column <- rlang::ensym(selected_column)
  plot_data <- data %>%
    mutate(
      selected = !!selected_column,
      plot_class = case_when(!selected ~ "Other", padj >= padj_cutoff ~ "Other", TRUE ~ metabolic_class)
    )

  n_selected <- sum(plot_data$selected)
  n_up <- sum(plot_data$selected & plot_data$padj < padj_cutoff & plot_data$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(plot_data$selected & plot_data$padj < padj_cutoff & plot_data$log2FoldChange < 0, na.rm = TRUE)

  cat(
    "\n============================================\n", output_paths$label, "\n",
    "============================================\n",
    "Selected genes : ", n_selected, "\n",
    "Up             : ", n_up, "\n",
    "Down           : ", n_down, "\n",
    sep = ""
  )

  p <- ggplot(plot_data, aes(x = log2FoldChange, y = y_plot)) +
    geom_point(
      data = plot_data %>% filter(plot_class == "Other", !above_y_cap),
      colour = "grey78", alpha = 0.55, size = 0.7
    ) +
    geom_point(
      data = plot_data %>% filter(plot_class != "Other", !above_y_cap),
      aes(colour = plot_class), alpha = 0.95, size = 1.5
    ) +
    geom_point(
      data = plot_data %>% filter(plot_class == "Other", above_y_cap),
      colour = "grey55", shape = 17, size = 2.2, show.legend = FALSE
    ) +
    geom_point(
      data = plot_data %>% filter(plot_class != "Other", above_y_cap),
      aes(colour = plot_class), shape = 17, size = 3, show.legend = FALSE
    ) +
    geom_point(
      data = apelin_points, shape = 21, fill = "white", colour = "black",
      stroke = 0.9, size = 2.8, show.legend = FALSE
    ) +
    scale_colour_manual(
      values = c("Lipid" = "darkorange2", "Amino acid" = "dodgerblue3", "Glucose" = "forestgreen"),
      breaks = c("Lipid", "Amino acid", "Glucose"), name = "Metabolic pathway"
    ) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", linewidth = 0.5) +
    scale_x_continuous(
      limits = c(-x_limit, x_limit), breaks = scales::breaks_pretty(n = 5),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(limits = c(0, y_cap), breaks = y_breaks, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 14) +
    theme(
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(colour = "black", size = 11, face = "bold"),
      axis.line = element_line(linewidth = 0.9, colour = "black"),
      axis.ticks = element_line(linewidth = 0.9, colour = "black"),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 10, b = 35, l = 10)
    ) +
    labs(x = expression(log[2]~"(Fold Change)"), y = expression(-log[10]~"(FDR)"))

  ggsave(output_paths$png, p, width = 6.4, height = 4.8, dpi = 600, bg = "white")
  ggsave(output_paths$pdf, p, width = 6.4, height = 4.8, bg = "white")
  ggsave(output_paths$svg, p, width = 6.4, height = 4.8, bg = "white")

  output_table <- plot_data %>%
    filter(selected) %>%
    arrange(padj) %>%
    dplyr::select(ENSEMBL, SYMBOL, baseMean, log2FoldChange, stat, pvalue, padj, lipid, amino_acid, glucose, metabolic_class)

  write_csv(output_table, output_paths$genes)
  p
}

# Build category sets
plot_metabolism <- plot_metabolism %>%
  mutate(
    set_lipid = lipid,
    set_lipid_amino = lipid | amino_acid,
    set_full_metabolism = lipid | amino_acid | glucose
  )

p_lipid <- make_metabolic_volcano(plot_metabolism, set_lipid, list(
  label = "03_metabolic_volcano_lipid",
  png = snakemake@output[["lipid_png"]], pdf = snakemake@output[["lipid_pdf"]],
  svg = snakemake@output[["lipid_svg"]], genes = snakemake@output[["lipid_genes"]]
))

p_lipid_amino <- make_metabolic_volcano(plot_metabolism, set_lipid_amino, list(
  label = "04_metabolic_volcano_lipid_aminoacid",
  png = snakemake@output[["lipid_amino_png"]], pdf = snakemake@output[["lipid_amino_pdf"]],
  svg = snakemake@output[["lipid_amino_svg"]], genes = snakemake@output[["lipid_amino_genes"]]
))

p_full <- make_metabolic_volcano(plot_metabolism, set_full_metabolism, list(
  label = "05_metabolic_volcano_lipid_aminoacid_glucose",
  png = snakemake@output[["full_png"]], pdf = snakemake@output[["full_pdf"]],
  svg = snakemake@output[["full_svg"]], genes = snakemake@output[["full_genes"]]
))

print(p_lipid)
print(p_lipid_amino)
print(p_full)

category_summary <- tibble(
  category = c("Lipid", "Amino acid", "Glucose", "Lipid + amino acid", "Lipid + amino acid + glucose"),
  n_genes = c(
    sum(plot_metabolism$lipid), sum(plot_metabolism$amino_acid), sum(plot_metabolism$glucose),
    sum(plot_metabolism$set_lipid_amino), sum(plot_metabolism$set_full_metabolism)
  )
)

print(category_summary)
write_csv(category_summary, snakemake@output[["summary"]])

cat(
  "\n============================================\n",
  "METABOLIC VOLCANO PLOTS COMPLETE\n",
  "============================================\n",
  "Colours:\n",
  "  Lipid       = orange\n",
  "  Amino acid  = blue\n",
  "  Glucose     = green\n\n",
  "Outputs:\n",
  "03 = Lipids\n",
  "04 = Lipids + amino acids\n",
  "05 = Lipids + amino acids + glucose\n",
  sep = ""
)
