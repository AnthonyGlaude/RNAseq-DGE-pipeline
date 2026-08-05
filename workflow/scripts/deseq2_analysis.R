suppressPackageStartupMessages({
  library(readr)
  library(DESeq2)
})

# — Snakemake inputs —
counts_file <- snakemake@input[["counts"]]
design_file <- snakemake@input[["metadata"]]

filter_thr <- as.integer(
  snakemake@params[["filter_count_threshold"]]
)

out_dir <- snakemake@params[["outdir"]]

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. Lecture de la matrice featureCounts fusionnée
# ============================================================

counts_df <- read_tsv(
  counts_file,
  col_types = cols()
)

# Colonnes produites par featureCounts qui ne sont pas des comptes
annotation_cols <- c(
  "Geneid",
  "Chr",
  "Start",
  "End",
  "Strand",
  "Length"
)

sample_cols <- setdiff(
  colnames(counts_df),
  annotation_cols
)

# Matrice numérique : gènes x échantillons
count_matrix <- as.matrix(
  counts_df[, sample_cols, drop = FALSE]
)

rownames(count_matrix) <- counts_df$Geneid

storage.mode(count_matrix) <- "integer"

# ============================================================
# 2. Lecture du design expérimental
# ============================================================

design_df <- read_tsv(
  design_file,
  col_types = cols()
)

# Vérifie que les IDs correspondent aux colonnes featureCounts
if (!all(design_df$sample_id %in% colnames(count_matrix))) {
  missing_samples <- setdiff(
    design_df$sample_id,
    colnames(count_matrix)
  )

  stop(
    "Échantillons absents de la matrice de comptes : ",
    paste(missing_samples, collapse = ", ")
  )
}

# Remet les colonnes dans le même ordre que le design
count_matrix <- count_matrix[
  ,
  design_df$sample_id,
  drop = FALSE
]

# Références :
# condition : chow
# genotype  : WT
design_df$condition <- factor(
  design_df$condition,
  levels = c("chow", "patho")
)

design_df$genotype <- factor(
  design_df$genotype,
  levels = c("WT", "KO")
)

sampleTable <- data.frame(
  condition = design_df$condition,
  genotype  = design_df$genotype,
  row.names = design_df$sample_id
)

# Vérification indispensable
stopifnot(
  identical(
    rownames(sampleTable),
    colnames(count_matrix)
  )
)

# ============================================================
# 3. Construction de l'objet DESeq2
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = sampleTable,
  design    = ~ condition * genotype
)

# ============================================================
# 4. Filtrage des faibles comptes
# ============================================================

keep <- rowSums(counts(dds)) >= filter_thr

dds <- dds[keep, ]

# ============================================================
# 5. Exécution de DESeq2
# ============================================================

dds <- DESeq(dds)

print(resultsNames(dds))

# ============================================================
# 6. Résultats principaux
# ============================================================

# Effet patho chez les WT :
# HFpEF_WT vs Control_WT
res_patho_WT <- results(
  dds,
  name = "condition_patho_vs_chow"
)

write.csv(
  as.data.frame(res_patho_WT),
  file = file.path(
    out_dir,
    "HFpEF_WT-Control_WT_DESeq2_gene.csv"
  ),
  quote = FALSE,
  row.names = TRUE
)

# Effet KO sous chow :
# Control_KO vs Control_WT
res_KO_chow <- results(
  dds,
  name = "genotype_KO_vs_WT"
)

write.csv(
  as.data.frame(res_KO_chow),
  file = file.path(
    out_dir,
    "Control_KO-Control_WT_DESeq2_gene.csv"
  ),
  quote = FALSE,
  row.names = TRUE
)

# Interaction :
# Est-ce que l'effet patho diffère entre KO et WT?
interaction_name <- grep(
  "condition.*genotype|genotype.*condition",
  resultsNames(dds),
  value = TRUE
)

if (length(interaction_name) != 1) {
  stop(
    "Impossible d'identifier automatiquement le coefficient ",
    "d'interaction. Coefficients disponibles : ",
    paste(resultsNames(dds), collapse = ", ")
  )
}

res_interaction <- results(
  dds,
  name = interaction_name
)

write.csv(
  as.data.frame(res_interaction),
  file = file.path(
    out_dir,
    "condition_genotype_interaction_DESeq2_gene.csv"
  ),
  quote = FALSE,
  row.names = TRUE
)

# Effet patho chez les KO :
# HFpEF_KO vs Control_KO
res_patho_KO <- results(
  dds,
  contrast = list(
    c(
      "condition_patho_vs_chow",
      interaction_name
    )
  )
)

write.csv(
  as.data.frame(res_patho_KO),
  file = file.path(
    out_dir,
    "HFpEF_KO-Control_KO_DESeq2_gene.csv"
  ),
  quote = FALSE,
  row.names = TRUE
)

# Effet KO dans la condition pathologique :
# HFpEF_KO vs HFpEF_WT
res_KO_patho <- results(
  dds,
  contrast = list(
    c(
      "genotype_KO_vs_WT",
      interaction_name
    )
  )
)

write.csv(
  as.data.frame(res_KO_patho),
  file = file.path(
    out_dir,
    "HFpEF_KO-HFpEF_WT_DESeq2_gene.csv"
  ),
  quote = FALSE,
  row.names = TRUE
)

saveRDS(
  dds,
  file = file.path(
    out_dir,
    "dds.rds"
  )
)
