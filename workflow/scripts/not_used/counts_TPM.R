suppressPackageStartupMessages({
  library(DESeq2)
  library(openxlsx)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

# ============================================================
# 1. Chemins
# ============================================================

dds_file <- paste0(
  "/home/glaudea/scratch/glaudea/test_souris_rnaseq/",
  "RNA_seq-analysis/workflow/results/deseq2/dds.rds"
)

featurecounts_file <- paste0(
  "/home/glaudea/scratch/glaudea/test_souris_rnaseq/",
  "RNA_seq-analysis/workflow/results/featurecounts/",
  "Control_WT_counts.txt"
)

output_file <- paste0(
  "/home/glaudea/scratch/glaudea/test_souris_rnaseq/",
  "RNA_seq-analysis/workflow/results/",
  "gene_expression_mean_counts_TPM.xlsx"
)

# ============================================================
# 2. Charger DDS
# ============================================================

dds <- readRDS(dds_file)

cat("Nombre de gènes :", nrow(dds), "\n")
cat("Nombre d'échantillons :", ncol(dds), "\n")

# ============================================================
# 3. Raw counts
# ============================================================

raw_counts <- counts(
  dds,
  normalized = FALSE
)

# ============================================================
# 4. Counts normalisés DESeq2
# ============================================================

normalized_counts <- counts(
  dds,
  normalized = TRUE
)

# ============================================================
# 5. Metadata / groupes
# ============================================================

metadata <- as.data.frame(colData(dds))

cat("\nColonnes metadata :\n")
print(colnames(metadata))

if ("group" %in% colnames(metadata)) {

  groups <- as.character(metadata$group)

} else if (
  "condition" %in% colnames(metadata) &&
  "genotype" %in% colnames(metadata)
) {

  groups <- paste(
    metadata$condition,
    metadata$genotype,
    sep = "_"
  )

} else {

  stop("Impossible de déterminer les groupes dans colData(dds).")
}

names(groups) <- colnames(dds)

cat("\nGroupes détectés :\n")
print(groups)

cat("\nNombre d'échantillons par groupe :\n")
print(table(groups))

# ============================================================
# 6. Ordre réel des groupes dans ton DDS
# ============================================================

group_order <- c(
  "chow_WT",
  "chow_KO",
  "patho_WT",
  "patho_KO"
)

# ============================================================
# 7. Fonction moyenne par groupe
# ============================================================

mean_by_group <- function(mat, groups, group_order) {

  result <- sapply(
    group_order,
    function(g) {

      samples <- names(groups)[groups == g]

      if (length(samples) == 0) {
        stop("Aucun échantillon trouvé pour le groupe : ", g)
      }

      rowMeans(
        mat[, samples, drop = FALSE],
        na.rm = TRUE
      )
    }
  )

  rownames(result) <- rownames(mat)
  colnames(result) <- group_order

  return(result)
}

# ============================================================
# 8. Moyennes raw counts
# ============================================================

raw_mean <- mean_by_group(
  raw_counts,
  groups,
  group_order
)

# ============================================================
# 9. Moyennes normalized counts
# ============================================================

normalized_mean <- mean_by_group(
  normalized_counts,
  groups,
  group_order
)

# ============================================================
# 10. Lire featureCounts
# ============================================================

fc <- read.delim(
  featurecounts_file,
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\nColonnes featureCounts :\n")
print(colnames(fc))

if (!"Geneid" %in% colnames(fc)) {
  stop("Colonne Geneid absente du fichier featureCounts.")
}

if (!"Length" %in% colnames(fc)) {
  stop("Colonne Length absente du fichier featureCounts.")
}

# ============================================================
# 11. Nettoyer IDs Ensembl
# ============================================================

clean_ensembl <- function(x) {
  sub("\\..*$", "", x)
}

dds_ids_clean <- clean_ensembl(
  rownames(raw_counts)
)

fc_ids_clean <- clean_ensembl(
  fc$Geneid
)

# ============================================================
# 12. Longueur des gènes
# ============================================================

gene_length_lookup <- fc$Length
names(gene_length_lookup) <- fc_ids_clean

gene_lengths <- gene_length_lookup[
  match(
    dds_ids_clean,
    names(gene_length_lookup)
  )
]

names(gene_lengths) <- rownames(raw_counts)

cat(
  "\nLongueurs trouvées :",
  sum(!is.na(gene_lengths)),
  "/",
  length(gene_lengths),
  "\n"
)

# ============================================================
# 13. Retirer les gènes sans longueur pour TPM
# ============================================================

valid_length <- (
  !is.na(gene_lengths) &
  gene_lengths > 0
)

cat(
  "Gènes utilisables pour TPM :",
  sum(valid_length),
  "\n"
)

# ============================================================
# 14. Calcul TPM par échantillon
# ============================================================

tpm <- matrix(
  NA_real_,
  nrow = nrow(raw_counts),
  ncol = ncol(raw_counts),
  dimnames = dimnames(raw_counts)
)

length_kb <- gene_lengths[valid_length] / 1000

rpk <- sweep(
  raw_counts[valid_length, , drop = FALSE],
  MARGIN = 1,
  STATS = length_kb,
  FUN = "/"
)

scaling_factor <- colSums(
  rpk,
  na.rm = TRUE
) / 1e6

tpm_valid <- sweep(
  rpk,
  MARGIN = 2,
  STATS = scaling_factor,
  FUN = "/"
)

tpm[valid_length, ] <- tpm_valid

# ============================================================
# 15. Vérification TPM
# ============================================================

cat("\nSomme TPM par échantillon :\n")
print(
  round(
    colSums(tpm, na.rm = TRUE)
  )
)

# Doit être environ 1 000 000 pour chaque échantillon.

# ============================================================
# 16. Moyenne TPM par groupe
# ============================================================

tpm_mean <- mean_by_group(
  tpm,
  groups,
  group_order
)

# ============================================================
# 17. Annotation SYMBOL
# ============================================================

symbols <- mapIds(
  org.Mm.eg.db,
  keys = dds_ids_clean,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

symbols <- unname(symbols)

# ============================================================
# 18. Fonction tableaux Excel
# ============================================================

make_output_table <- function(mat) {

  data.frame(
    ENSEMBL_ID = rownames(mat),
    SYMBOL = symbols,

    Control_WT = mat[, "chow_WT"],
    Control_KO = mat[, "chow_KO"],
    HFpEF_WT = mat[, "patho_WT"],
    HFpEF_KO = mat[, "patho_KO"],

    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 19. Tables finales
# ============================================================

raw_df <- make_output_table(
  raw_mean
)

normalized_df <- make_output_table(
  normalized_mean
)

tpm_df <- make_output_table(
  tpm_mean
)

# ============================================================
# 20. Créer Excel
# ============================================================

wb <- createWorkbook()

addWorksheet(
  wb,
  "Raw_counts_mean"
)

writeData(
  wb,
  "Raw_counts_mean",
  raw_df,
  withFilter = TRUE
)

addWorksheet(
  wb,
  "Normalized_counts_mean"
)

writeData(
  wb,
  "Normalized_counts_mean",
  normalized_df,
  withFilter = TRUE
)

addWorksheet(
  wb,
  "TPM_mean"
)

writeData(
  wb,
  "TPM_mean",
  tpm_df,
  withFilter = TRUE
)

# ============================================================
# 21. Mise en forme
# ============================================================

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center"
)

for (sheet in names(wb)) {

  freezePane(
    wb,
    sheet,
    firstRow = TRUE,
    firstCol = TRUE
  )

  addStyle(
    wb,
    sheet = sheet,
    style = header_style,
    rows = 1,
    cols = 1:6,
    gridExpand = TRUE
  )

  setColWidths(
    wb,
    sheet,
    cols = 1,
    widths = 24
  )

  setColWidths(
    wb,
    sheet,
    cols = 2,
    widths = 16
  )

  setColWidths(
    wb,
    sheet,
    cols = 3:6,
    widths = 16
  )
}

# ============================================================
# 22. Sauvegarder
# ============================================================

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

saveWorkbook(
  wb,
  output_file,
  overwrite = TRUE
)

cat("\n============================================\n")
cat("Excel créé :\n")
cat(output_file, "\n")
cat("============================================\n")

cat("\nColonnes finales :\n")
cat(
  "ENSEMBL_ID | SYMBOL | Control_WT | Control_KO | HFpEF_WT | HFpEF_KO\n"
)
