## merge_counts.R
## Fusionne les fichiers de comptage featureCounts (un par condition, generes
## par la rule feature_counts_by_condition) en une seule matrice de comptage,
## prete pour DESeq2.
##
## Chaque fichier d'entree contient :
##   - 1re ligne : ligne de commande featureCounts (commence par '#')
##   - 2e ligne : header (Geneid, Chr, Start, End, Strand, Length, <bam paths>)
## Les colonnes de metadonnees (Chr, Start, End, Strand, Length) sont
## identiques dans les 4 fichiers puisqu'ils partagent la meme annotation GTF ;
## on ne les garde qu'une seule fois dans la sortie finale.

args <- commandArgs(trailingOnly = TRUE)
input_files <- args[-length(args)]
output_file <- args[length(args)]

## Colonnes de metadonnees generees par featureCounts, en plus de Geneid
meta_cols <- c("Chr", "Start", "End", "Strand", "Length")

merged <- NULL

for (f in input_files) {

  df <- read.table(
    f,
    header = TRUE,
    sep = "\t",
    skip = 1,          # saute la ligne de commande featureCounts
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  ## Renomme les colonnes d'echantillons : featureCounts utilise le chemin
  ## complet du BAM comme nom de colonne (ex: results/STAR/sample1/Aligned...)
  ## On ne garde que l'id de l'echantillon, extrait du chemin.
  sample_cols <- setdiff(colnames(df), c("Geneid", meta_cols))
  new_names <- basename(dirname(sample_cols))
  colnames(df)[colnames(df) %in% sample_cols] <- new_names

  if (is.null(merged)) {
    ## Premier fichier : on garde Geneid + metadonnees + comptages
    merged <- df
  } else {
    ## Fichiers suivants : on ne garde que Geneid + comptages, on verifie
    ## que les metadonnees matchent (meme GTF, meme ordre de genes attendu)
    stopifnot(identical(df$Geneid, merged$Geneid))

    df_counts <- df[, c("Geneid", new_names), drop = FALSE]
    merged <- merge(merged, df_counts, by = "Geneid", sort = FALSE)
  }
}

## Verification finale : pas de colonne d'echantillon dupliquee
sample_columns_final <- setdiff(colnames(merged), c("Geneid", meta_cols))
if (any(duplicated(sample_columns_final))) {
  stop("Des noms d'echantillons sont dupliques apres fusion : ",
       paste(sample_columns_final[duplicated(sample_columns_final)], collapse = ", "))
}

write.table(
  merged,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Matrice fusionnee :", nrow(merged), "genes x",
    length(sample_columns_final), "echantillons\n")
cat("Echantillons :", paste(sample_columns_final, collapse = ", "), "\n")
