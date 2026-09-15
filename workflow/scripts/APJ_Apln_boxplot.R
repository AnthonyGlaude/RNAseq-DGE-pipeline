suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(openxlsx)
  library(ragg)
  library(svglite)
})

# Inputs and outputs
dds_file <- snakemake@input[["dds"]]
stats_file <- snakemake@input[["stats"]]
out_dir <- dirname(snakemake@output[["png"]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Liberation Sans"
genes_of_interest <- data.frame(
  gene_id = c("ENSMUSG00000044338", "ENSMUSG00000037010"),
  gene_name = c("Aplnr", "Apln"),
  stringsAsFactors = FALSE
)

# Load the WT-only DDS and the corresponding DESeq2 results
dds <- readRDS(dds_file)
stats <- read_csv(stats_file, show_col_types = FALSE)

cat("\nSamples in DDS:\n")
print(colnames(dds))
cat("\nConditions:\n")
print(table(colData(dds)$condition))

if (ncol(dds) != 8) warning("The DDS contains ", ncol(dds), " samples instead of the 8 expected WT samples.")
if (!"condition" %in% colnames(colData(dds))) stop("The DDS does not contain a 'condition' column.")

required_stats <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
missing_stats <- setdiff(required_stats, names(stats))
if (length(missing_stats) > 0) stop("Missing DESeq2 columns: ", paste(missing_stats, collapse = ", "))

gene_candidates <- c("gene", "gene_id", "ensembl_gene_id", "ENSEMBL", "...1")
gene_col <- gene_candidates[gene_candidates %in% names(stats)][1]
if (is.na(gene_col)) stop("Could not identify the gene ID column in: ", paste(names(stats), collapse = ", "))
stats <- stats %>% mutate(ENSEMBL = sub("\\..*$", "", as.character(.data[[gene_col]])))

get_sig_label <- function(padj) {
  if (is.na(padj)) return("ns")
  if (padj < 0.0001) return("****")
  if (padj < 0.001) return("***")
  if (padj < 0.01) return("**")
  if (padj < 0.05) return("*")
  "ns"
}

standardize_condition <- function(x) {
  case_when(
    tolower(x) %in% c("chow", "control", "control_wt") ~ "Chow",
    tolower(x) %in% c("patho", "hfpef", "hfpef_wt") ~ "HFpEF",
    TRUE ~ NA_character_
  )
}

plot_list <- list()
stats_list <- list()

for (i in seq_len(nrow(genes_of_interest))) {
  gene_id <- genes_of_interest$gene_id[i]
  gene_name <- genes_of_interest$gene_name[i]
  cat("\nProcessing: ", gene_name, " (", gene_id, ")\n", sep = "")

  if (!gene_id %in% rownames(dds)) stop("Gene ", gene_id, " is absent from the DDS.")
  gene_stats <- stats %>% filter(ENSEMBL == gene_id)
  if (nrow(gene_stats) != 1) stop("Expected one DESeq2 row for ", gene_id, ", found ", nrow(gene_stats), ".")

  padj <- gene_stats$padj[[1]]
  sig_label <- get_sig_label(padj)
  df <- plotCounts(dds, gene = gene_id, intgroup = "condition", returnData = TRUE)
  df$sample <- rownames(df)
  df$condition <- standardize_condition(as.character(df$condition))
  if (anyNA(df$condition)) stop("Unrecognized condition values in DDS: ", paste(unique(colData(dds)$condition), collapse = ", "))

  df$condition <- factor(df$condition, levels = c("Chow", "HFpEF"))
  df$log2_count <- log2(df$count + 1)
  df$gene_id <- gene_id
  df$gene_name <- gene_name
  df$padj <- padj
  df$sig_label <- sig_label
  plot_list[[gene_name]] <- df

  stats_list[[gene_name]] <- data.frame(
    gene = gene_name,
    gene_id = gene_id,
    comparison = "HFpEF WT vs Chow WT",
    baseMean = gene_stats$baseMean[[1]],
    log2FoldChange = gene_stats$log2FoldChange[[1]],
    lfcSE = gene_stats$lfcSE[[1]],
    stat = gene_stats$stat[[1]],
    pvalue = gene_stats$pvalue[[1]],
    padj = padj,
    sig_label = sig_label,
    row.names = NULL
  )

  cat("FDR = ", signif(padj, 4), " -> ", sig_label, "\n", sep = "")
}

plot_df <- bind_rows(plot_list)
stats_df <- bind_rows(stats_list)
df_export <- plot_df %>%
  transmute(gene = gene_name, gene_id, sample, count, condition = as.character(condition), log2_count)

cat("\nValues used for the figure:\n")
print(df_export)

# Export normalized counts and DESeq2 statistics
wb <- createWorkbook()
addWorksheet(wb, "Normalized_counts")
writeData(wb, "Normalized_counts", df_export)
setColWidths(wb, "Normalized_counts", cols = seq_len(ncol(df_export)), widths = "auto")
addWorksheet(wb, "DESeq2_statistics")
writeData(wb, "DESeq2_statistics", stats_df)
setColWidths(wb, "DESeq2_statistics", cols = seq_len(ncol(stats_df)), widths = "auto")
saveWorkbook(wb, snakemake@output[["xlsx"]], overwrite = TRUE)

# Plot positions
plot_df$gene_name <- factor(plot_df$gene_name, levels = c("Aplnr", "Apln"))
plot_df$group <- interaction(plot_df$gene_name, plot_df$condition, sep = "_")
plot_df$x_pos <- case_when(
  plot_df$gene_name == "Aplnr" & plot_df$condition == "Chow" ~ 1.00,
  plot_df$gene_name == "Aplnr" & plot_df$condition == "HFpEF" ~ 1.45,
  plot_df$gene_name == "Apln" & plot_df$condition == "Chow" ~ 2.15,
  plot_df$gene_name == "Apln" & plot_df$condition == "HFpEF" ~ 2.60
)

aplnr_sig <- stats_df$sig_label[stats_df$gene == "Aplnr"]
apln_sig <- stats_df$sig_label[stats_df$gene == "Apln"]
y_max <- max(plot_df$log2_count, na.rm = TRUE)
y_min <- min(plot_df$log2_count, na.rm = TRUE)
y_range <- max(y_max - y_min, 1)
bracket_y <- y_max + 0.08 * y_range
tip_length <- 0.04 * y_range
text_y <- y_max + 0.13 * y_range
plot_ymin <- floor((y_min - 0.08 * y_range) * 2) / 2
plot_ymax <- ceiling((text_y + 0.08 * y_range) * 2) / 2
y_breaks <- seq(plot_ymin, plot_ymax, by = 0.5)

p <- ggplot(plot_df, aes(x = x_pos, y = log2_count, fill = condition)) +
  geom_boxplot(aes(group = group, colour = condition), width = 0.22, outlier.shape = NA, linewidth = 0.9) +
  geom_jitter(aes(colour = condition, group = group), width = 0.025, height = 0, size = 3, alpha = 1) +
  annotate("segment", x = 1.00, xend = 1.45, y = bracket_y, yend = bracket_y, linewidth = 0.8) +
  annotate("segment", x = 1.00, xend = 1.00, y = bracket_y, yend = bracket_y - tip_length, linewidth = 0.8) +
  annotate("segment", x = 1.45, xend = 1.45, y = bracket_y, yend = bracket_y - tip_length, linewidth = 0.8) +
  annotate("text", x = 1.225, y = text_y, label = aplnr_sig, family = font_family, fontface = "bold", size = 5) +
  annotate("segment", x = 2.15, xend = 2.60, y = bracket_y, yend = bracket_y, linewidth = 0.8) +
  annotate("segment", x = 2.15, xend = 2.15, y = bracket_y, yend = bracket_y - tip_length, linewidth = 0.8) +
  annotate("segment", x = 2.60, xend = 2.60, y = bracket_y, yend = bracket_y - tip_length, linewidth = 0.8) +
  annotate("text", x = 2.375, y = text_y, label = apln_sig, family = font_family, fontface = "bold", size = 5) +
  scale_fill_manual(values = c("Chow" = "#A0A0A440", "HFpEF" = "#FF808040")) +
  scale_colour_manual(values = c("Chow" = "#A0A0A4", "HFpEF" = "#FF8080")) +
  scale_x_continuous(
    limits = c(0.70, 2.90), breaks = c(1.00, 1.45, 2.15, 2.60),
    labels = c("Chow", "HFpEF", "Chow", "HFpEF"), expand = c(0, 0)
  ) +
  scale_y_continuous(breaks = y_breaks, expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = expression(log[2] * "(DESeq2 normalized counts + 1)")) +
  coord_cartesian(ylim = c(plot_ymin, plot_ymax), clip = "off") +
  theme_classic(base_size = 14, base_family = font_family) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold", size = 13, colour = "black"),
    axis.text.x = element_text(face = "bold", size = 12, colour = "black"),
    axis.text.y = element_text(face = "bold", size = 11, colour = "black"),
    axis.line = element_line(linewidth = 0.9, colour = "black"),
    axis.ticks = element_line(linewidth = 0.9, colour = "black"),
    plot.margin = margin(t = 20, r = 10, b = 35, l = 10)
  ) +
  annotate("text", x = 1.225, y = -Inf, label = "Aplnr", vjust = 4, fontface = "bold", family = font_family, size = 4.5) +
  annotate("text", x = 2.375, y = -Inf, label = "Apln", vjust = 4, fontface = "bold", family = font_family, size = 4.5)

ggsave(snakemake@output[["png"]], p, device = ragg::agg_png, width = 6.4, height = 4.8, units = "in", res = 600, bg = "white")
ggsave(snakemake@output[["pdf"]], p, device = grDevices::cairo_pdf, width = 6.4, height = 4.8, units = "in", bg = "white")
ggsave(snakemake@output[["svg"]], p, device = svglite::svglite, width = 6.4, height = 4.8, units = "in", bg = "white")

cat(
  "\n============================================\n",
  "Analysis complete\n",
  "============================================\n",
  "PNG  : ", snakemake@output[["png"]], "\n",
  "PDF  : ", snakemake@output[["pdf"]], "\n",
  "SVG  : ", snakemake@output[["svg"]], "\n",
  "XLSX : ", snakemake@output[["xlsx"]], "\n",
  sep = ""
)
