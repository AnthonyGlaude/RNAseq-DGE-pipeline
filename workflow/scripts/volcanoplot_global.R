#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# Inputs and outputs
stats_file <- snakemake@input[["stats"]]
volcano_dir <- dirname(snakemake@output[["fdr_png"]])
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters
padj_cutoff <- 0.05
fc_cutoff <- 1.2
lfc_cutoff <- log2(fc_cutoff)
y_cap <- 80
y_cap_position <- 78
y_breaks <- c(0, 20, 40, 60, 80)

cat(
  "\n============================================\n",
  "VOLCANO PARAMETERS\n",
  "============================================\n",
  "FDR cutoff       :", padj_cutoff, "\n",
  "FC cutoff        :", fc_cutoff, "\n",
  "log2FC cutoff    :", round(lfc_cutoff, 3), "\n",
  "Visual Y cap     :", y_cap, "\n"
)

# Read and validate DESeq2 results
df <- read_csv(stats_file, show_col_types = FALSE)
required_columns <- c("log2FoldChange", "padj")
missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

# Prepare plotting data
res_plot <- df %>%
  filter(!is.na(log2FoldChange), !is.na(padj), padj > 0) %>%
  mutate(
    y_real = -log10(padj),
    above_y_cap = y_real > y_cap,
    y_plot = if_else(above_y_cap, y_cap_position, y_real)
  )

cat("\nGenes included in volcano plots:", nrow(res_plot), "\n")
cat("Genes above Y-axis cap:", sum(res_plot$above_y_cap), "\n")

cat(
  "\n============================================\n",
  "DATA RANGE\n",
  "============================================\n",
  "Minimum log2FC :", round(min(res_plot$log2FoldChange, na.rm = TRUE), 3), "\n",
  "Maximum log2FC :", round(max(res_plot$log2FoldChange, na.rm = TRUE), 3), "\n",
  "Maximum -log10(FDR) :", round(max(res_plot$y_real, na.rm = TRUE), 3), "\n"
)

# Volcano 01: FDR only
plot_01 <- res_plot %>%
  mutate(volcano_class = case_when(
    padj < padj_cutoff & log2FoldChange > 0 ~ "Up",
    padj < padj_cutoff & log2FoldChange < 0 ~ "Down",
    TRUE ~ "NS"
  ))

n_up_01 <- sum(plot_01$volcano_class == "Up")
n_down_01 <- sum(plot_01$volcano_class == "Down")
n_sig_01 <- n_up_01 + n_down_01

cat(
  "\n============================================\n",
  "VOLCANO 01 - FDR ONLY\n",
  "============================================\n",
  "Up   :", n_up_01, "\n",
  "Down :", n_down_01, "\n",
  "Total significant :", n_sig_01, "\n"
)

# Volcano 02: FDR < 0.05 and FC >= 1.2
plot_02 <- res_plot %>%
  mutate(volcano_class = case_when(
    padj < padj_cutoff & log2FoldChange >= lfc_cutoff ~ "Up",
    padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "Down",
    TRUE ~ "NS"
  ))

n_up_02 <- sum(plot_02$volcano_class == "Up")
n_down_02 <- sum(plot_02$volcano_class == "Down")
n_sig_02 <- n_up_02 + n_down_02

cat(
  "\n============================================\n",
  "VOLCANO 02 - FDR + FC ", fc_cutoff, "\n",
  "============================================\n",
  "Up   :", n_up_02, "\n",
  "Down :", n_down_02, "\n",
  "Total significant :", n_sig_02, "\n"
)

make_volcano <- function(plot_data, n_up, n_down, add_fc_cutoff = FALSE) {
  max_abs_x <- max(abs(plot_data$log2FoldChange), na.rm = TRUE)
  x_limit <- ceiling(max_abs_x * 1.05)
  if (x_limit < 1) x_limit <- 1

  p <- ggplot(plot_data, aes(x = log2FoldChange, y = y_plot)) +
    geom_point(
      data = plot_data %>% filter(volcano_class == "NS", !above_y_cap),
      colour = "grey75", alpha = 0.60, size = 0.7
    ) +
    geom_point(
      data = plot_data %>% filter(volcano_class != "NS", !above_y_cap),
      aes(colour = volcano_class), alpha = 0.90, size = 0.8
    ) +
    geom_point(
      data = plot_data %>% filter(volcano_class == "NS", above_y_cap),
      colour = "grey55", shape = 17, size = 2.5, show.legend = FALSE
    ) +
    geom_point(
      data = plot_data %>% filter(volcano_class != "NS", above_y_cap),
      aes(colour = volcano_class), shape = 17, size = 3, show.legend = FALSE
    ) +
    scale_colour_manual(values = c("Up" = "red2", "Down" = "dodgerblue3")) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", linewidth = 0.5) +
    annotate(
      "text", x = -0.92 * x_limit, y = 75, label = paste0("Down: ", n_down),
      hjust = 0, colour = "dodgerblue3", fontface = "bold", size = 3.8
    ) +
    annotate(
      "text", x = 0.92 * x_limit, y = 75, label = paste0("Up: ", n_up),
      hjust = 1, colour = "red2", fontface = "bold", size = 3.8
    )

  if (add_fc_cutoff) {
    p <- p + geom_vline(
      xintercept = c(-lfc_cutoff, lfc_cutoff),
      linetype = "dashed", linewidth = 0.5
    )
  }

  p +
    scale_x_continuous(
      limits = c(-x_limit, x_limit),
      breaks = scales::breaks_pretty(n = 5),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(limits = c(0, y_cap), breaks = y_breaks, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 11) +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(colour = "black", size = 9),
      legend.position = "none",
      plot.margin = margin(t = 8, r = 8, b = 5, l = 8)
    ) +
    labs(x = expression(log[2]~"(Fold Change)"), y = expression(-log[10]~"(FDR)"))
}

p_01 <- make_volcano(plot_01, n_up_01, n_down_01, add_fc_cutoff = FALSE)
p_02 <- make_volcano(plot_02, n_up_02, n_down_02, add_fc_cutoff = TRUE)

print(p_01)
print(p_02)

# Save plots using the paths declared in the Snakemake rule
ggsave(snakemake@output[["fdr_png"]], p_01, width = 6, height = 5, dpi = 600, bg = "white")
ggsave(snakemake@output[["fdr_pdf"]], p_01, width = 6, height = 5, bg = "white")
ggsave(snakemake@output[["fdr_svg"]], p_01, width = 6, height = 5, bg = "white")
ggsave(snakemake@output[["fc_png"]], p_02, width = 6, height = 5, dpi = 600, bg = "white")
ggsave(snakemake@output[["fc_pdf"]], p_02, width = 6, height = 5, bg = "white")
ggsave(snakemake@output[["fc_svg"]], p_02, width = 6, height = 5, bg = "white")

# Export gene counts
summary_counts <- tibble(
  volcano = c("01_FDR_only", "02_FDR_FC1.2"),
  FDR = c(padj_cutoff, padj_cutoff),
  FC_cutoff = c(NA, fc_cutoff),
  log2FC_cutoff = c(NA, lfc_cutoff),
  Up = c(n_up_01, n_up_02),
  Down = c(n_down_01, n_down_02),
  Total_significant = c(n_sig_01, n_sig_02)
)

print(summary_counts)
write_csv(summary_counts, snakemake@output[["counts"]])

cat(
  "\n============================================\n",
  "VOLCANO PLOTS COMPLETE\n",
  "============================================\n",
  "Output directory:", volcano_dir, "\n"
)
