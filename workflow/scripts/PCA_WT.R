suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(ragg)
  library(svglite)
})

dds_file <- snakemake@input[["dds"]]
out_dir <- dirname(snakemake@output[["png"]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ntop <- 1000
dds <- readRDS(dds_file)

message("Samples:")
print(colnames(dds))
message("\nMetadata:")
print(as.data.frame(colData(dds)))

# Assign groups from sample names
condition_pca <- case_when(
  grepl("^Control_WT", colnames(dds)) ~ "Control",
  grepl("^HFpEF_WT", colnames(dds)) ~ "HFpEF",
  TRUE ~ NA_character_
)

if (anyNA(condition_pca)) {
  stop("Samples not recognized as Control_WT or HFpEF_WT: ", paste(colnames(dds)[is.na(condition_pca)], collapse = ", "))
}

condition_pca <- factor(condition_pca, levels = c("Control", "HFpEF"))
if (ncol(dds) != 8) warning("The DDS contains ", ncol(dds), " samples instead of the 8 expected WT samples.")
if (any(table(condition_pca) < 2)) stop("PCA requires at least two samples in each condition.")

message("\nPCA groups:")
print(table(condition_pca))

# Variance-stabilizing transformation and selection of the 1000 most variable genes
vsd <- vst(dds, blind = TRUE)
vst_matrix <- assay(vsd)
gene_variances <- apply(vst_matrix, 1, var)
n_selected <- min(ntop, length(gene_variances))
selected_genes <- order(gene_variances, decreasing = TRUE)[seq_len(n_selected)]
vst_top <- vst_matrix[selected_genes, , drop = FALSE]

# PCA
pca <- prcomp(t(vst_top), center = TRUE, scale. = FALSE)
percent_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_df <- as.data.frame(pca$x) %>%
  mutate(
    sample = rownames(pca$x),
    condition = condition_pca[match(sample, colnames(dds))],
    sample_label = sub("^.*_", "", sample)
  ) %>%
  relocate(sample, sample_label, condition)

variance_df <- tibble(
  component = paste0("PC", seq_along(percent_var)),
  variance_percent = percent_var,
  cumulative_percent = cumsum(percent_var)
)

selected_df <- tibble(
  gene_id = rownames(vst_matrix)[selected_genes],
  variance = gene_variances[selected_genes]
)

write_csv(pca_df, snakemake@output[["coordinates"]])
write_csv(variance_df, snakemake@output[["variance"]])
write_csv(selected_df, snakemake@output[["selected_genes"]])

message("\nPCA coordinates:")
print(pca_df %>% select(sample, condition, PC1, PC2))

condition_colors <- c("Control" = "#0072B2", "HFpEF" = "#D55E00")
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample_label)) +
  geom_point(size = 4.5, alpha = 0.9) +
  geom_text_repel(size = 4, show.legend = FALSE) +
  scale_color_manual(values = condition_colors) +
  labs(
    title = "PCA - Control WT vs HFpEF WT",
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)"),
    color = "Condition"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_pca)
ggsave(snakemake@output[["png"]], p_pca, device = ragg::agg_png, width = 8, height = 6.5, units = "in", dpi = 300, bg = "white")
ggsave(snakemake@output[["pdf"]], p_pca, device = grDevices::cairo_pdf, width = 8, height = 6.5, units = "in", bg = "white")
ggsave(snakemake@output[["svg"]], p_pca, device = svglite::svglite, width = 8, height = 6.5, units = "in", bg = "white")

message("\nPCA saved: ", snakemake@output[["png"]])
message("Genes used: ", n_selected)
message("PC1 = ", round(percent_var[1], 1), "% ; PC2 = ", round(percent_var[2], 1), "%")
