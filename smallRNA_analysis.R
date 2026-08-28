# ==============================================================================
# Project: Small RNA-seq analysis of sperm in the equine epididymis
#
# Description:
#   Preprocessing and filtering of small RNA-seq count data.
#   This script performs sample matching, expression filtering,
#   and downstream analyses.
#
# Author: Youness EL AMMANY
# Affiliation: CR2TI
# Created: 2026-08-28
#
# Input:
#   - count_matrix.tsv
#   - metadata.tsv
#
#
# R version: 4.5.3
# Dependencies:
#   CRAN:
#     tidyverse 2.0.0
#     dplyr 1.2.0
#     tidyr 1.3.2
#     tibble 3.3.1
#     readr 2.2.0
#     ggplot2 4.0.3
#     ggrepel 0.9.8
#     VennDiagram 1.8.2
#     ggforce 0.5.0
#     pheatmap 1.0.13
#     pdftools 3.9.0
#     knitr 1.51
#
#   Bioconductor:
#     DESeq2 1.50.2
#     SummarizedExperiment 1.40.0
#     BiocManager 1.30.27
#
# Repository:
#   https://github.com/YounsElmny/smallRNA-pipeline
#
# ==============================================================================


## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = TRUE,
  warning = FALSE,
  fig.align = "center",
  fig.width = 12,
  fig.height = 8,
  out.width = "100%"
)


## ----helper-include-pdf-as-large-png, include=FALSE----------------------------------------------------------------------------------
if (!requireNamespace("pdftools", quietly = TRUE)) {
  install.packages("pdftools", repos = "https://cloud.r-project.org")
}

include_pdf_as_large_png <- function(figs, dpi = 220, width = "100%") {
  
  figs <- figs[file.exists(figs)]
  
  if (length(figs) == 0) {
    cat("<p><em>Aucune figure trouvée pour ce bloc.</em></p>")
    return(invisible(NULL))
  }
  
  png_files <- gsub("\\.pdf$", ".png", figs, ignore.case = TRUE)
  
  for (i in seq_along(figs)) {
    pdf_file <- figs[i]
    png_file <- png_files[i]
    
    if (!file.exists(png_file) || file.info(png_file)$mtime < file.info(pdf_file)$mtime) {
      message("Conversion PDF vers PNG : ", pdf_file)
      
      invisible(utils::capture.output(
        suppressMessages(suppressWarnings(
          pdftools::pdf_convert(
            pdf = pdf_file,
            format = "png",
            dpi = dpi,
            filenames = png_file
          )
        )),
        type = "output"
      ))
    }
  }
  
  png_files <- png_files[file.exists(png_files)]
  
  for (png in png_files) {
    img_uri <- knitr::image_uri(png)
    cat(
      '<div class="figure-block">',
      '<img src="', img_uri, '" style="width:', width, ';">',
      '</div>',
      sep = ""
    )
  }
  
  invisible(NULL)
}


## ----1-packages-colors-paths, eval=TRUE, message=FALSE, warning=FALSE, results='hide'------------------------------------------------
############################################################
# 1. Packages, colors, paths
############################################################

## Explicitly set a CRAN mirror to allow package installation
## automatically when knitting the RMarkdown document.
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = getOption("repos")[["CRAN"]])
}

cran_pkgs <- c(
  "tidyverse", "readr", "tibble", "ggrepel", "ggplot2",
  "VennDiagram", "ggforce"
)
bioc_pkgs <- c("DESeq2", "SummarizedExperiment")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = getOption("repos")[["CRAN"]])
  }
}

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}

library(tidyverse)
library(ggrepel)
library(DESeq2)
library(SummarizedExperiment)

## Colors
cols_region <- c(
  Head = "#F1D21A",
  Body = "#F28E1C",
  Tail = "#D90D0D",
  `N/A` = "#D4A017"
)

cols_tissue <- c(
  Epididymis = "#4F81BD",
  Spermatozoa = "#9DC3E6",
  Testis = "#2F5597"
)

cols_treatment <- c(
  `Anti-GnRH` = "#BDBDBD",
  None = "#D9B44A"
)

cols_biotype <- c(
  "miR" = "#566D7E",
  "tRF" = "#78909C",
  "rRF" = "#9AAEB8",
  "snoRNA" = "#A9897D",
  "piR" = "#8A817C",
  "ncRNA" = "#AAA39D",
  "cdna" = "#C4B8AE",
  "Unknown" = "#D9D9D9"
)

cols_volcano <- c(
  "Not_DE" = "#C7C7C7",
  "Higher_in_Head" = "#617C93",
  "Higher_in_Body" = "#617C93",
  "Higher_in_Tail" = "#A77768"
)

## Paths
setwd("WorkDirectory")

counts_file <- "Count_Matrix.tsv"
meta_file   <- "QNS31219_metadata.txt"

outdir <- "Figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## General thresholds
padj_cut <- 0.05
lfc_cut  <- 1
topN     <- 15


## ----2-read-counts-metadata, eval=TRUE-----------------------------------------------------------------------------------------------
############################################################
# 2. Read counts + metadata
############################################################

counts <- read.delim(
  counts_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = "",
  fileEncoding = "Windows-1252"
)

meta <- read.delim(
  meta_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = "",
  dec = ",",
  fileEncoding = "Windows-1252"
)

## Remove any empty columns
counts <- counts[, colnames(counts) != "", drop = FALSE]


## ----3-metadata-cleaning-sample-columns, eval=TRUE-----------------------------------------------------------------------------------
############################################################
# 3. Metadata cleaning + sample columns
############################################################

required_count_cols <- c("sequence")
required_meta_cols <- c(
  "nb",
  "Unique Sample name",
  "Origin in Epididymis",
  "Tissue type",
  "Treatment"
)

missing_count_cols <- setdiff(required_count_cols, colnames(counts))
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))

if (length(missing_count_cols) > 0) {
  stop("Colonnes manquantes dans counts : ", paste(missing_count_cols, collapse = ", "))
}

if (length(missing_meta_cols) > 0) {
  stop("Colonnes manquantes dans metadata : ", paste(missing_meta_cols, collapse = ", "))
}

## Sample columns A01-A42
sample_cols <- grep("^A[0-9]{2}$", colnames(counts), value = TRUE)

if (length(sample_cols) == 0) {
  stop("Aucune colonne sample de type A01, A02, ..., A42 trouvée.")
}

message("Nombre de colonnes samples trouvées : ", length(sample_cols))
print(sample_cols)

## Raw count matrix
counts_only <- counts[, sample_cols, drop = FALSE]

counts_only <- apply(counts_only, 2, function(x) {
  x <- gsub('"', "", x)
  x <- trimws(x)
  as.integer(x)
})

counts_only <- as.matrix(counts_only)
rownames(counts_only) <- make.unique(counts$sequence)

if (anyNA(counts_only)) {
  warning("Il y a des NA après conversion en entier dans counts_only.")
}

## Match A01-A42 with metadata
meta$A_sample <- sprintf("A%02d", as.integer(meta$nb))

missing_in_counts <- setdiff(meta$A_sample, colnames(counts_only))
missing_in_meta   <- setdiff(colnames(counts_only), meta$A_sample)

if (length(missing_in_counts) > 0) {
  stop("Samples présents dans metadata mais absents des counts : ",
       paste(missing_in_counts, collapse = ", "))
}

if (length(missing_in_meta) > 0) {
  stop("Colonnes counts absentes du metadata : ",
       paste(missing_in_meta, collapse = ", "))
}

## Metadata cleaning
meta <- meta %>%
  dplyr::mutate(
    `Tissue type` = trimws(as.character(`Tissue type`)),
    `Origin in Epididymis` = trimws(as.character(`Origin in Epididymis`)),
    Treatment = trimws(as.character(Treatment))
  )

meta$`Origin in Epididymis`[
  is.na(meta$`Origin in Epididymis`) |
    meta$`Origin in Epididymis` == "" |
    meta$`Origin in Epididymis` == "N/A"
] <- "N/A"

meta$`Origin in Epididymis` <- factor(
  meta$`Origin in Epididymis`,
  levels = c("Head", "Body", "Tail", "N/A")
)

meta$`Tissue type` <- factor(
  meta$`Tissue type`,
  levels = c("Epididymis", "Spermatozoa", "Testis")
)

meta$Treatment <- factor(
  meta$Treatment,
  levels = c("Anti-GnRH", "None")
)

## Reorder metadata to match counts_only
idx <- match(colnames(counts_only), meta$A_sample)

if (any(is.na(idx))) {
  stop("Samples counts absents de meta$A_sample : ",
       paste(colnames(counts_only)[is.na(idx)], collapse = ", "))
}

meta <- meta[idx, , drop = FALSE]
rownames(meta) <- meta$A_sample

stopifnot(all(rownames(meta) == colnames(counts_only)))


## ----4-expression-filtering, eval=TRUE----
############################################################
# 4. Expression filtering
############################################################

message("---- Filtrage expression ----")

min_count <- 10

## Groups = major tissue types
group_samples <- split(meta$A_sample, meta$`Tissue type`)
group_samples <- lapply(group_samples, function(x) intersect(x, colnames(counts_only)))

group_sizes <- lengths(group_samples)
min_samp_by_group <- ceiling(group_sizes * 0.5)

filter_summary <- tibble::tibble(
  tissue_group = names(group_samples),
  n_samples = group_sizes,
  min_samples_required = min_samp_by_group
)

print(filter_summary)

write.table(
  filter_summary,
  file = file.path(outdir, "filter_tissue_groups_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

keep_by_group <- sapply(names(group_samples), function(g) {
  samples_g <- group_samples[[g]]
  min_samp_g <- min_samp_by_group[[g]]
  rowSums(counts_only[, samples_g, drop = FALSE] >= min_count) >= min_samp_g
})

keep_by_group <- as.matrix(keep_by_group)
keep <- rowSums(keep_by_group) >= 1

counts_filt <- counts_only[keep, , drop = FALSE]
counts_filtered_out <- counts[keep, , drop = FALSE]

message("Séquences gardées après filtre : ", nrow(counts_filt), " / ", nrow(counts_only))
message("Séquences retirées après filtre : ", sum(!keep))

write.table(
  counts_filtered_out,
  file = file.path(outdir, "42_samples_matrix_count.filtered_before_PCA.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## ----5-utility-functions, eval=TRUE--------------------------------------------------------------------------------------------------
############################################################
# 5. Utility functions
############################################################

empty_to_na <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "na", "n/a", "-", ".", "NULL", "null")] <- NA
  x
}

is_filled <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) & !(x %in% c("", "NA", "N/A", "na", "n/a", "-", ".", "NULL", "null"))
}

pick_one_random_annotation <- function(x) {
  if (is.na(x)) return(NA_character_)
  x <- trimws(as.character(x))
  if (x == "") return("")

  annot <- unlist(strsplit(x, split = ","))
  annot <- trimws(annot)
  annot <- annot[annot != ""]

  if (length(annot) == 0) return("")
  sample(annot, size = 1)
}

clean_feature_text <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- "unannotated"
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}


clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

safe_vst <- function(dds_obj, label = "") {
  tryCatch(
    {
      message("VST avec DESeq2::vst() ", label)
      DESeq2::vst(dds_obj, blind = TRUE)
    },
    error = function(e) {
      message("DESeq2::vst() a échoué ", label)
      message("Erreur : ", conditionMessage(e))
      message("Utilisation de varianceStabilizingTransformation().")
      DESeq2::varianceStabilizingTransformation(dds_obj, blind = TRUE)
    }
  )
}

prepare_spz_metadata <- function(regions = c("Head", "Body", "Tail")) {
  meta_spz <- meta %>%
    dplyr::filter(
      `Tissue type` == "Spermatozoa",
      `Origin in Epididymis` %in% regions
    )

  meta_spz <- as.data.frame(meta_spz)
  rownames(meta_spz) <- meta_spz$A_sample

  meta_spz$Origin_in_Epididymis <- as.factor(meta_spz$`Origin in Epididymis`)
  meta_spz$Origin_in_Epididymis <- droplevels(meta_spz$Origin_in_Epididymis)

  return(meta_spz)
}


## ----6-random-single-annotation-per-annotation-cell, eval=TRUE------
############################################################
# 6. Random single annotation per annotation cell
############################################################

message("---- Simplification annotations multiples ----")

set.seed(42)

cols_random_annot <- c(
  "matures",
  "Rfam_snoRNA_eca",
  "Rfam_rRNA_eca",
  "pirbase_eca_v30",
  "gtrnadb_mature_tRNA_eca",
  "ENSEMBL_eca_ncrna",
  "ENSEMBL_eca_cdna"
)

cols_random_annot <- intersect(cols_random_annot, colnames(counts_filtered_out))

message("Colonnes simplifiées : ", paste(cols_random_annot, collapse = ", "))

counts_filtered_out[cols_random_annot] <- lapply(
  counts_filtered_out[cols_random_annot],
  function(col) vapply(col, pick_one_random_annotation, character(1))
)

write.table(
  counts_filtered_out,
  file = file.path(outdir, "42_samples_matrix_count.filtered_single_random_annotation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## ----7-priority-annotation-gene-name-priority-biotype-priority, eval=TRUE-----
############################################################
# 7. Priority annotation: gene_name_priority + biotype_priority
############################################################

message("---- Annotation prioritaire ----")

required_annot_cols <- c(
  "matures",
  "gtrnadb_mature_tRNA_eca",
  "Rfam_rRNA_eca",
  "Rfam_snoRNA_eca",
  "pirbase_eca_v30",
  "ENSEMBL_eca_ncrna",
  "ENSEMBL_eca_cdna"
)

missing_annot_cols <- setdiff(required_annot_cols, colnames(counts_filtered_out))

if (length(missing_annot_cols) > 0) {
  stop("Colonnes d'annotation manquantes : ",
       paste(missing_annot_cols, collapse = ", "))
}

## Clean the columns used for priority annotation
for (col in required_annot_cols) {
  counts_filtered_out[[col]] <- empty_to_na(counts_filtered_out[[col]])
}

## Create the priority annotation
## Priority: miR > tRF > rRF > snoRNA > piR > ncRNA > cdna
counts_filtered_out <- counts_filtered_out %>%
  dplyr::mutate(
    gene_name_priority = dplyr::case_when(
      !is.na(matures) ~ matures,
      !is.na(gtrnadb_mature_tRNA_eca) ~ gtrnadb_mature_tRNA_eca,
      !is.na(Rfam_rRNA_eca) ~ Rfam_rRNA_eca,
      !is.na(Rfam_snoRNA_eca) ~ Rfam_snoRNA_eca,
      !is.na(pirbase_eca_v30) ~ pirbase_eca_v30,
      !is.na(ENSEMBL_eca_ncrna) ~ ENSEMBL_eca_ncrna,
      !is.na(ENSEMBL_eca_cdna) ~ ENSEMBL_eca_cdna,
      TRUE ~ NA_character_
    ),
    biotype_priority = dplyr::case_when(
      !is.na(matures) ~ "miR",
      !is.na(gtrnadb_mature_tRNA_eca) ~ "tRF",
      !is.na(Rfam_rRNA_eca) ~ "rRF",
      !is.na(Rfam_snoRNA_eca) ~ "snoRNA",
      !is.na(pirbase_eca_v30) ~ "piR",
      !is.na(ENSEMBL_eca_ncrna) ~ "ncRNA",
      !is.na(ENSEMBL_eca_cdna) ~ "cdna",
      TRUE ~ NA_character_
    )
  )

priority_summary <- counts_filtered_out %>%
  dplyr::group_by(biotype_priority) %>%
  dplyr::summarise(n_sequences = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_sequences))

print(priority_summary)

write.table(
  priority_summary,
  file = file.path(outdir, "annotation_priority_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  counts_filtered_out,
  file = file.path(outdir, "42_samples_matrix_count.filtered_with_priority_annotation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## Unannotated diagnostic
still_unannotated_check <- counts_filtered_out %>%
  dplyr::filter(is.na(gene_name_priority) | gene_name_priority == "") %>%
  dplyr::select(
    sequence,
    matures,
    Rfam_snoRNA_eca,
    Rfam_rRNA_eca,
    pirbase_eca_v30,
    gtrnadb_mature_tRNA_eca,
    ENSEMBL_eca_ncrna,
    ENSEMBL_eca_cdna,
    dplyr::everything()
  )

message("Nombre de lignes encore sans gene_name_priority : ", nrow(still_unannotated_check))

write.table(
  still_unannotated_check,
  file = file.path(outdir, "check_still_unannotated_after_priority.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## ----8-read-level-feature-id-counts-de, eval=TRUE------------------------------------------------------------------------------------
############################################################
# 8. Read-level feature_id + counts_DE
############################################################

message("---- Création feature_id read-level ----")

required_feature_cols <- c("sequence", "gene_name_priority")
missing_feature_cols <- setdiff(required_feature_cols, colnames(counts_filtered_out))

if (length(missing_feature_cols) > 0) {
  stop("Colonnes manquantes pour créer feature_id : ",
       paste(missing_feature_cols, collapse = ", "))
}

counts_filtered_out <- counts_filtered_out %>%
  dplyr::mutate(
    sequence_clean = clean_feature_text(sequence),
    gene_name_clean = clean_feature_text(gene_name_priority),
    feature_id = paste(sequence_clean, gene_name_clean, sep = "_")
  )

if (nrow(counts_filtered_out) != length(unique(counts_filtered_out$feature_id))) {
  warning("Certains feature_id sont dupliqués. Ajout d'un suffixe unique avec make.unique().")
  counts_filtered_out$feature_id <- make.unique(counts_filtered_out$feature_id, sep = "_dup")
}

stopifnot(length(unique(counts_filtered_out$feature_id)) == nrow(counts_filtered_out))

counts_DE <- counts_filtered_out[, sample_cols, drop = FALSE]

counts_DE <- apply(counts_DE, 2, function(x) {
  x <- gsub('"', "", x)
  x <- trimws(x)
  as.integer(x)
})

counts_DE <- as.matrix(counts_DE)
rownames(counts_DE) <- counts_filtered_out$feature_id

stopifnot(all(colnames(counts_DE) == rownames(meta)))

write.table(
  counts_filtered_out,
  file = file.path(outdir, "42_samples_matrix_count.filtered_with_feature_id.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

counts_DE_out <- as.data.frame(counts_DE) %>%
  tibble::rownames_to_column("feature_id")

write.table(
  counts_DE_out,
  file = file.path(outdir, "counts_matrix_for_DESeq2_feature_id.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Dimensions counts_DE read-level : ", paste(dim(counts_DE), collapse = " x "))


## ----9-gene-level-matrix-counts-gene-de-annot-gene, eval=TRUE------------------------------------------------------------------------
############################################################
# 9. Gene-level matrix: counts_gene_DE + annot_gene
############################################################

message("---- Création matrice gene-level ----")

counts_gene_input <- counts_filtered_out %>%
  dplyr::filter(
    !is.na(gene_name_priority),
    gene_name_priority != "",
    gene_name_priority != "unannotated"
  ) %>%
  dplyr::mutate(
    gene_id = clean_gene_id(gene_name_priority)
  )

message("Nombre de reads/séquences annotés utilisés : ", nrow(counts_gene_input))
message("Nombre de gènes uniques avant nettoyage : ", length(unique(counts_gene_input$gene_name_priority)))
message("Nombre de gene_id uniques après nettoyage : ", length(unique(counts_gene_input$gene_id)))

annot_gene <- counts_gene_input %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(
    gene_name_priority = dplyr::first(gene_name_priority),
    biotype_priority = dplyr::first(biotype_priority),
    n_reads_merged = dplyr::n(),
    .groups = "drop"
  )

counts_gene_df <- counts_gene_input %>%
  dplyr::select(gene_id, dplyr::all_of(sample_cols)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(sample_cols),
      ~ sum(as.integer(.x), na.rm = TRUE)
    ),
    .groups = "drop"
  )

counts_gene_DE <- counts_gene_df %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

############## Check ##############
#counts_gene_input %>%
#  dplyr::filter(gene_id == "CM009148_1_112302588_112302054") %>%
#  dplyr::select(
#    sequence,
#    gene_name_priority,
#    biotype_priority,
#    dplyr::all_of(sample_cols)
#  )
#########################################

storage.mode(counts_gene_DE) <- "integer"

## Check count conservation
total_before_merge <- sum(as.matrix(counts_gene_input[, sample_cols]), na.rm = TRUE)
total_after_merge  <- sum(counts_gene_DE, na.rm = TRUE)

message("Total counts avant merge : ", total_before_merge)
message("Total counts après merge : ", total_after_merge)

if (total_before_merge == total_after_merge) {
  message("OK : le total des counts est conservé après merge.")
} else {
  warning("Attention : le total des counts diffère après merge.")
}

write.table(
  counts_gene_df,
  file = file.path(outdir, "counts_matrix_gene_level_merged_reads.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  annot_gene,
  file = file.path(outdir, "annotation_gene_level_merged_reads.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

top_merged_genes <- annot_gene %>%
  dplyr::arrange(dplyr::desc(n_reads_merged)) %>%
  dplyr::slice_head(n = 20)

write.table(
  top_merged_genes,
  file = file.path(outdir, "check_top20_genes_most_reads_merged.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Dimensions counts_gene_DE : ", paste(dim(counts_gene_DE), collapse = " x "))


## ----10-qc-biotype-proportions-histograms, eval=TRUE------
############################################################
# 10. QC / biotype proportions / histograms
############################################################

message("---- QC et proportions ----")

## Library size
lib_sizes <- colSums(counts_only)

qc_lib <- tibble::tibble(
  A_sample = names(lib_sizes),
  lib_size = as.numeric(lib_sizes)
) %>%
  dplyr::left_join(meta %>% as.data.frame(), by = "A_sample")

pdf(file.path(outdir, "QC_lib_size_42_samples_matrix.pdf"),
    width = 10, height = 9, useDingbats = FALSE)

ggplot(qc_lib, aes(x = reorder(`Unique Sample name`, lib_size), y = lib_size)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Taille de librairie - matrice 42 samples",
    x = NULL,
    y = "Total reads"
  ) +
  theme_bw()

invisible(dev.off())

## Distribution of the number of distinct reads per gene
gene_read_distribution <- counts_filtered_out %>%
  dplyr::filter(!is.na(gene_name_priority), gene_name_priority != "") %>%
  dplyr::group_by(gene_name_priority, biotype_priority) %>%
  dplyr::summarise(
    n_different_reads = dplyr::n_distinct(sequence),
    total_counts = sum(rowSums(dplyr::across(dplyr::all_of(sample_cols))), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_different_reads))

write.table(
  gene_read_distribution,
  file = file.path(outdir, "gene_read_distribution.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

gene_read_distribution_plot <- gene_read_distribution %>%
  dplyr::mutate(
    n_reads_category = dplyr::case_when(
      n_different_reads > 10 ~ ">10",
      TRUE ~ as.character(n_different_reads)
    ),
    n_reads_category = factor(n_reads_category, levels = c(as.character(1:10), ">10"))
  )

pdf(file.path(outdir, "Histogram_number_of_different_reads_per_gene_grouped_gt10.pdf"),
    width = 9, height = 7, useDingbats = FALSE)

ggplot(gene_read_distribution_plot, aes(x = n_reads_category)) +
  geom_bar(color = "black", fill = "grey70") +
  labs(
    title = "Distribution du nombre de reads différents par gène",
    x = "Nombre de reads/séquences différentes par gène",
    y = "Nombre de gènes"
  ) +
  theme_bw()

invisible(dev.off())

## Biotype proportions before merging: number of distinct reads/sequences
biotype_prop_before_merge <- counts_filtered_out %>%
  dplyr::filter(!is.na(biotype_priority), biotype_priority != "") %>%
  dplyr::group_by(biotype_priority) %>%
  dplyr::summarise(n_reads_different = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(
    proportion = n_reads_different / sum(n_reads_different),
    percent = 100 * proportion
  ) %>%
  dplyr::arrange(dplyr::desc(n_reads_different))

write.table(
  biotype_prop_before_merge,
  file = file.path(outdir, "biotype_proportions_before_merge_by_number_of_reads.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf(file.path(outdir, "Biotype_proportions_before_merge_by_number_of_reads.pdf"),
    width = 9, height = 7, useDingbats = FALSE)

ggplot(
  biotype_prop_before_merge,
  aes(x = reorder(biotype_priority, -percent), y = percent)
) +
  geom_col(color = "black", fill = "grey70") +
  geom_text(aes(label = paste0(round(percent, 1), "%")),
            vjust = -0.3, size = 3.5) +
  labs(
    title = "Proportions des biotypes avant merge",
    subtitle = "Niveau read/séquence : 1 ligne = 1 read différent",
    x = "Biotype prioritaire",
    y = "Pourcentage de reads/séquences différentes"
  ) +
  theme_bw()

invisible(dev.off())

## Biotype proportions based on total counts
counts_biotype_input <- counts_filtered_out %>%
  dplyr::filter(!is.na(biotype_priority), biotype_priority != "") %>%
  dplyr::mutate(
    total_reads_all_samples = rowSums(
      dplyr::across(dplyr::all_of(sample_cols)),
      na.rm = TRUE
    )
  )

biotype_prop_counts <- counts_biotype_input %>%
  dplyr::group_by(biotype_priority) %>%
  dplyr::summarise(
    total_reads = sum(total_reads_all_samples, na.rm = TRUE),
    n_different_reads = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    proportion = total_reads / sum(total_reads),
    percent = 100 * proportion
  ) %>%
  dplyr::arrange(dplyr::desc(total_reads))

write.table(
  biotype_prop_counts,
  file = file.path(outdir, "biotype_proportions_after_priority_by_total_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf(file.path(outdir, "Biotype_proportions_after_priority_by_total_counts.pdf"),
    width = 9, height = 7, useDingbats = FALSE)

ggplot(
  biotype_prop_counts,
  aes(x = reorder(biotype_priority, -percent), y = percent)
) +
  geom_col(color = "black", fill = "grey70") +
  geom_text(aes(label = paste0(round(percent, 1), "%")),
            vjust = -0.3, size = 3.5) +
  labs(
    title = "Proportions des biotypes selon le nombre total de reads",
    subtitle = "Après gene priority : chaque read/séquence appartient à un seul biotype",
    x = "Biotype prioritaire",
    y = "Pourcentage du total des reads comptés"
  ) +
  theme_bw()

invisible(dev.off())

## Distribution of reads by biotype and sample category
## Representation harmonized with the other analyses: pie charts + stacked barplots

## Define spermatozoa samples, also used below for the cDNA check ----
sperm_samples <- meta %>%
  dplyr::filter(`Tissue type` == "Spermatozoa") %>%
  dplyr::pull(A_sample) %>%
  intersect(sample_cols)

## Clean / prepare sample categories ----
meta_biotype <- meta %>%
  dplyr::distinct(A_sample, .keep_all = TRUE) %>%
  dplyr::mutate(
    A_sample = trimws(as.character(A_sample)),
    `Tissue type` = trimws(as.character(`Tissue type`)),
    `Origin in Epididymis` = trimws(as.character(`Origin in Epididymis`)),
    `Tissue type` = dplyr::case_when(
      is.na(`Tissue type`) | `Tissue type` == "" ~ "Unknown",
      TRUE ~ `Tissue type`
    ),
    `Origin in Epididymis` = dplyr::case_when(
      is.na(`Origin in Epididymis`) |
        `Origin in Epididymis` == "" |
        `Origin in Epididymis` == "N/A" ~ "NA",
      TRUE ~ `Origin in Epididymis`
    ),
    sample_category = paste(`Tissue type`, `Origin in Epididymis`, sep = "_")
  )

## Calculate reads by biotype and sample ----
## First aggregate by biotype to avoid creating an excessively large table.
biotype_by_sample <- counts_filtered_out %>%
  dplyr::mutate(
    Biotype = dplyr::case_when(
      is.na(biotype_priority) | biotype_priority == "" ~ "Unknown",
      TRUE ~ as.character(biotype_priority)
    )
  ) %>%
  dplyr::group_by(Biotype) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(sample_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(sample_cols),
    names_to = "A_sample",
    values_to = "reads_biotype"
  ) %>%
  dplyr::mutate(
    A_sample = trimws(as.character(A_sample)),
    reads_biotype = as.numeric(reads_biotype)
  ) %>%
  dplyr::left_join(
    meta_biotype %>%
      dplyr::select(A_sample, `Tissue type`, `Origin in Epididymis`, Treatment, sample_category),
    by = "A_sample"
  ) %>%
  dplyr::group_by(A_sample) %>%
  dplyr::mutate(
    total_reads_sample = sum(reads_biotype, na.rm = TRUE),
    percent_reads_sample = 100 * reads_biotype / total_reads_sample
  ) %>%
  dplyr::ungroup()

write.table(
  biotype_by_sample,
  file = file.path(outdir, "Biotype_read_percentages_by_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## Calculate reads by biotype and category ----
biotype_by_category <- biotype_by_sample %>%
  dplyr::group_by(sample_category, `Tissue type`, `Origin in Epididymis`, Biotype) %>%
  dplyr::summarise(
    reads_biotype = sum(reads_biotype, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(sample_category) %>%
  dplyr::mutate(
    total_reads_category = sum(reads_biotype, na.rm = TRUE),
    percent_reads_category = 100 * reads_biotype / total_reads_category
  ) %>%
  dplyr::ungroup()

write.table(
  biotype_by_category,
  file = file.path(outdir, "Biotype_read_percentages_by_category.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## Version for figures: no grouping into Other ----
biotype_by_category_plot <- biotype_by_category %>%
  dplyr::mutate(Biotype_plot = Biotype)

## Faceted pie charts by category ----
pdf(
  file.path(outdir, "Piecharts_biotype_percent_reads_by_category.pdf"),
  width = 14,
  height = 10,
  useDingbats = FALSE
)

ggplot(
  biotype_by_category_plot,
  aes(
    x = "",
    y = percent_reads_category,
    fill = Biotype_plot
  )
) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ sample_category, ncol = 3) +
  scale_fill_manual(
    values = cols_biotype,
    drop = FALSE
  ) +
  labs(
    title = "Pourcentage de reads par biotype",
    subtitle = "Catégories = Tissue type + Origin in Epididymis",
    fill = "Biotype"
  ) +
  theme_void() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

invisible(dev.off())

## Stacked percentage barplot, easier to read for comparing categories ----
pdf(
  file.path(outdir, "Barplot_biotype_percent_reads_by_category.pdf"),
  width = 14,
  height = 8,
  useDingbats = FALSE
)

ggplot(
  biotype_by_category_plot,
  aes(
    x = sample_category,
    y = percent_reads_category,
    fill = Biotype_plot
  )
) +
  geom_col(color = "white") +
  coord_flip() +
  labs(
    title = "Pourcentage de reads par biotype et par catégorie d'échantillons",
    subtitle = "Tous les biotypes sont affichés séparément",
    x = NULL,
    y = "% de reads",
    fill = "Biotype"
  ) +
  scale_fill_manual(
    values = cols_biotype,
    drop = FALSE
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold")
  )

invisible(dev.off())

## Additional barplot: raw read counts by biotype and category ----
pdf(
  file.path(outdir, "Barplot_biotype_raw_reads_by_category.pdf"),
  width = 14,
  height = 8,
  useDingbats = FALSE
)

ggplot(
  biotype_by_category_plot,
  aes(
    x = sample_category,
    y = reads_biotype,
    fill = Biotype_plot
  )
) +
  geom_col(color = "white") +
  coord_flip() +
  labs(
    title = "Quantité brute de reads par biotype et par catégorie d'échantillons",
    subtitle = "Reads = somme des counts par biotype",
    x = NULL,
    y = "Nombre de reads",
    fill = "Biotype"
  ) +
  scale_fill_manual(
    values = cols_biotype,
    drop = FALSE
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold")
  )

invisible(dev.off())

message("Analyse de répartition des biotypes terminée.")
message("Fichiers générés :")
message("- Biotype_read_percentages_by_sample.tsv")
message("- Biotype_read_percentages_by_category.tsv")
message("- Piecharts_biotype_percent_reads_by_category.pdf")
message("- Barplot_biotype_percent_reads_by_category.pdf")
message("- Barplot_biotype_raw_reads_by_category.pdf")

## Number of distinct cDNA reads in sperm
cdna_reads_sperm <- counts_filtered_out %>%
  dplyr::filter(biotype_priority == "cdna") %>%
  dplyr::mutate(
    total_reads_sperm = rowSums(
      dplyr::across(dplyr::all_of(sperm_samples)),
      na.rm = TRUE
    )
  ) %>%
  dplyr::filter(total_reads_sperm > 0)

n_cdna_different_reads_sperm <- nrow(cdna_reads_sperm)

message("Nombre de reads différents cDNA présents dans le sperm : ",
        n_cdna_different_reads_sperm)

write.table(
  cdna_reads_sperm,
  file = file.path(outdir, "cdna_different_reads_in_sperm.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## ----figures-10-qc-biotype-proportions-histograms, echo=FALSE, message=FALSE, warning=FALSE, results="asis"--------------------------
figs <- c(
  file.path(outdir, "QC_lib_size_42_samples_matrix.pdf"),
  file.path(outdir, "Histogram_number_of_different_reads_per_gene_grouped_gt10.pdf"),
  file.path(outdir, "Biotype_proportions_before_merge_by_number_of_reads.pdf"),
  file.path(outdir, "Biotype_proportions_after_priority_by_total_counts.pdf"),
  file.path(outdir, "Piecharts_biotype_percent_reads_by_category.pdf"),
  file.path(outdir, "Barplot_biotype_percent_reads_by_category.pdf"),
  file.path(outdir, "Barplot_biotype_raw_reads_by_category.pdf")
)
include_pdf_as_large_png(figs)


## ----11-global-pca, eval=TRUE--------------------------------------------------------------------------------------------------------
############################################################
# 11. Global PCA
############################################################

message("---- PCA globales ----")

make_pca_plot <- function(vst_matrix, meta_df, filename, color_by, shape_by = "Origin in Epididymis", title_prefix = "PCA") {

  pca <- prcomp(t(vst_matrix), scale. = FALSE)
  percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

  pca_df <- as.data.frame(pca$x[, 1:2]) %>%
    tibble::rownames_to_column("A_sample") %>%
    dplyr::left_join(meta_df %>% as.data.frame(), by = "A_sample")

  p <- ggplot(
    pca_df,
    aes(
      x = PC1,
      y = PC2,
      color = .data[[color_by]],
      label = `Unique Sample name`
    )
  )

  if (is.null(shape_by)) {
    p <- p + geom_point(size = 4, alpha = 0.9)
  } else {
    p <- p + geom_point(aes(shape = .data[[shape_by]]), size = 4, alpha = 0.9)
  }

  p <- p +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 100) +
    labs(
      title = title_prefix,
      x = paste0("PC1 : ", percentVar[1], "%"),
      y = paste0("PC2 : ", percentVar[2], "%"),
      color = color_by
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )

  if (!is.null(shape_by)) {
    p <- p + labs(shape = shape_by)
  }

  if (color_by == "Origin in Epididymis") {
    p <- p + scale_color_manual(values = cols_region, drop = FALSE)
  }

  if (color_by == "Treatment") {
    p <- p + scale_color_manual(values = cols_treatment, drop = FALSE)
  }

  pdf(file.path(outdir, filename), width = 9, height = 7, useDingbats = FALSE)
  print(p)
  invisible(dev.off())

  return(pca_df)
}

## VST PCA for all samples
dds_pca_all <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(counts_filt),
  colData = meta,
  design = ~ 1
)

vsd_pca_all <- safe_vst(dds_pca_all, label = "PCA all samples")
vst_mat <- SummarizedExperiment::assay(vsd_pca_all)


pca_all_region <- make_pca_plot(
  vst_matrix = vst_mat,
  meta_df = meta,
  filename = "PCA_all_samples_color_Origin_in_Epididymis.pdf",
  color_by = "Origin in Epididymis",
  shape_by = "Tissue type",
  title_prefix = "PCA VST all samples"
)

pca_all_treatment <- make_pca_plot(
  vst_matrix = vst_mat,
  meta_df = meta,
  filename = "PCA_all_samples_color_Treatment.pdf",
  color_by = "Treatment",
  shape_by = "Tissue type",
  title_prefix = "PCA VST all samples"
)


## PCA spermatozoa only
meta_spz <- meta %>% dplyr::filter(`Tissue type` == "Spermatozoa")
counts_spz <- counts_filt[, meta_spz$A_sample, drop = FALSE]

dds_spz <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(counts_spz),
  colData = meta_spz,
  design = ~ 1
)

vsd_spz <- safe_vst(dds_spz, label = "PCA sperm")
vst_spz <- SummarizedExperiment::assay(vsd_spz)

pca_spz_region <- make_pca_plot(
  vst_matrix = vst_spz,
  meta_df = meta_spz,
  filename = "PCA_spermatozoa_only_color_Origin_in_Epididymis.pdf",
  color_by = "Origin in Epididymis",
  shape_by = NULL,
  title_prefix = "PCA VST spermatozoa"
)

write.table(
  pca_spz_region,
  file = file.path(outdir, "PCA_coordinates_spermatozoa_only.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## PCA epididymis only
meta_ep <- meta %>% dplyr::filter(`Tissue type` == "Epididymis")
counts_ep <- counts_filt[, meta_ep$A_sample, drop = FALSE]

dds_ep <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(counts_ep),
  colData = meta_ep,
  design = ~ 1
)

vsd_ep <- safe_vst(dds_ep, label = "PCA epididymis")
vst_ep <- SummarizedExperiment::assay(vsd_ep)

pca_ep_region <- make_pca_plot(
  vst_matrix = vst_ep,
  meta_df = meta_ep,
  filename = "PCA_epididymis_only_color_Origin_in_Epididymis.pdf",
  color_by = "Origin in Epididymis",
  shape_by = NULL,
  title_prefix = "PCA VST epididymis"
)

write.table(
  pca_ep_region,
  file = file.path(outdir, "PCA_coordinates_epididymis_only.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## ----figures-11-global-pca, echo=FALSE, message=FALSE, warning=FALSE, results="asis"-------------------------------------------------
figs <- c(
  file.path(outdir, "PCA_all_samples_color_Tissue_type.pdf"),
  file.path(outdir, "PCA_all_samples_color_Origin_in_Epididymis.pdf"),
  file.path(outdir, "PCA_all_samples_color_Treatment.pdf"),
  file.path(outdir, "PCA_spermatozoa_only_color_Origin_in_Epididymis.pdf"),
  file.path(outdir, "PCA_epididymis_only_color_Origin_in_Epididymis.pdf"),
  file.path(outdir, "PCA_log_counts_simple.pdf")
)
include_pdf_as_large_png(figs)


## ----12-pca-by-biotype-in-spermatozoa, eval=TRUE-------------------------------------------------------------------------------------
############################################################
# 12. PCA by biotype in spermatozoa
############################################################

message("---- PCA par biotype dans les spermatozoïdes ----")

pca_biotype_dir <- file.path(outdir, "PCA_by_biotype_SPZ")
dir.create(pca_biotype_dir, showWarnings = FALSE, recursive = TRUE)

meta_spz_pca <- meta %>%
  dplyr::filter(`Tissue type` == "Spermatozoa")

meta_spz_pca <- as.data.frame(meta_spz_pca)
rownames(meta_spz_pca) <- meta_spz_pca$A_sample

spz_samples <- intersect(rownames(meta_spz_pca), colnames(counts_DE))
meta_spz_pca <- meta_spz_pca[spz_samples, , drop = FALSE]

message("Nombre de samples SPZ pour PCA : ", length(spz_samples))
print(table(meta_spz_pca$`Origin in Epididymis`, useNA = "ifany"))

make_pca_one_biotype <- function(biotype_name) {

  message("PCA biotype : ", biotype_name)

  features_biotype <- counts_filtered_out %>%
    dplyr::filter(biotype_priority == biotype_name) %>%
    dplyr::pull(feature_id) %>%
    intersect(rownames(counts_DE))

  message("Nombre de features pour ", biotype_name, " : ", length(features_biotype))

  if (length(features_biotype) < 2) {
    warning("Biotype ", biotype_name, " ignoré : moins de 2 features.")
    return(NULL)
  }

  counts_bio_spz <- counts_DE[features_biotype, spz_samples, drop = FALSE]
  counts_bio_spz <- counts_bio_spz[rowSums(counts_bio_spz) > 0, , drop = FALSE]

  message("Features gardées après filtre >0 : ", nrow(counts_bio_spz))

  if (nrow(counts_bio_spz) < 2) {
    warning("Biotype ", biotype_name, " ignoré après filtre : moins de 2 features.")
    return(NULL)
  }

  dds_bio <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts_bio_spz),
    colData = meta_spz_pca,
    design = ~ 1
  )

  vsd_bio <- safe_vst(dds_bio, label = paste0("PCA biotype ", biotype_name))
  vst_bio <- SummarizedExperiment::assay(vsd_bio)

  pca_bio <- prcomp(t(vst_bio), scale. = FALSE)
  percentVar <- round(100 * (pca_bio$sdev^2 / sum(pca_bio$sdev^2)), 1)

  pca_bio_df <- as.data.frame(pca_bio$x[, 1:2]) %>%
    tibble::rownames_to_column("A_sample") %>%
    dplyr::left_join(meta_spz_pca %>% as.data.frame(), by = "A_sample")

  write.table(
    pca_bio_df,
    file = file.path(pca_biotype_dir, paste0("PCA_coordinates_SPZ_", biotype_name, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  pdf(file.path(pca_biotype_dir, paste0("PCA_SPZ_", biotype_name, ".pdf")),
      width = 9, height = 7, useDingbats = FALSE)

  p <- ggplot(
    pca_bio_df,
    aes(
      x = PC1,
      y = PC2,
      color = `Origin in Epididymis`,
      label = `Unique Sample name`
    )
  ) +
    geom_point(size = 4, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 100) +
    scale_color_manual(values = cols_region, drop = FALSE) +
    labs(
      title = paste0("PCA SPZ — biotype ", biotype_name),
      subtitle = paste0(nrow(counts_bio_spz), " features utilisées"),
      x = paste0("PC1 : ", percentVar[1], "%"),
      y = paste0("PC2 : ", percentVar[2], "%"),
      color = "Origin in Epididymis"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )

  print(p)
  invisible(dev.off())

  return(pca_bio_df)
}

biotypes_to_pca <- counts_filtered_out %>%
  dplyr::filter(!is.na(biotype_priority), biotype_priority != "") %>%
  dplyr::pull(biotype_priority) %>%
  unique()

print(biotypes_to_pca)

pca_by_biotype_results <- lapply(biotypes_to_pca, make_pca_one_biotype)
names(pca_by_biotype_results) <- biotypes_to_pca

message("PCA par biotype terminée.")


## ----figures-12-pca-by-biotype-in-spermatozoa, echo=FALSE, message=FALSE, warning=FALSE, results="asis"------------------------------
fig_dir <- file.path(outdir, "PCA_by_biotype_SPZ")
figs <- if (dir.exists(fig_dir)) list.files(fig_dir, pattern = "\\.pdf$", full.names = TRUE) else character(0)
include_pdf_as_large_png(figs)


## ----13-pca-by-biotype-in-epididymis, eval=TRUE--------------------------------------------------------------------------------------
############################################################
# 13. PCA by biotype in epididymis
############################################################

message("---- PCA par biotype dans l'épididyme ----")

pca_biotype_ep_dir <- file.path(outdir, "PCA_by_biotype_Epididymis")
dir.create(pca_biotype_ep_dir, showWarnings = FALSE, recursive = TRUE)

## Epididymis-only metadata: Head / Body / Tail
meta_ep_pca <- meta %>%
  dplyr::filter(
    `Tissue type` == "Epididymis",
    `Origin in Epididymis` %in% c("Head", "Body", "Tail")
  )

meta_ep_pca <- as.data.frame(meta_ep_pca)
rownames(meta_ep_pca) <- meta_ep_pca$A_sample

## Epididymis samples present in the read-level matrix
ep_samples <- intersect(rownames(meta_ep_pca), colnames(counts_DE))
meta_ep_pca <- meta_ep_pca[ep_samples, , drop = FALSE]

## Explicitly preserve the Head / Body / Tail order
meta_ep_pca$`Origin in Epididymis` <- factor(
  meta_ep_pca$`Origin in Epididymis`,
  levels = c("Head", "Body", "Tail")
)

message("Nombre de samples épididyme pour PCA : ", length(ep_samples))
print(table(meta_ep_pca$`Origin in Epididymis`, useNA = "ifany"))

if (length(ep_samples) < 3) {
  stop("Moins de 3 échantillons d'épididyme disponibles : PCA impossible.")
}

## PCA function for one biotype in the epididymis
make_pca_one_biotype_epididymis <- function(biotype_name) {

  message("PCA épididyme — biotype : ", biotype_name)

  ## Features belonging to the requested biotype
  features_biotype <- counts_filtered_out %>%
    dplyr::filter(biotype_priority == biotype_name) %>%
    dplyr::pull(feature_id) %>%
    intersect(rownames(counts_DE))

  message("Nombre de features pour ", biotype_name, " : ", length(features_biotype))

  if (length(features_biotype) < 2) {
    warning("Biotype ", biotype_name, " ignoré : moins de 2 features.")
    return(NULL)
  }

  ## Biotype matrix restricted to epididymis samples
  counts_bio_ep <- counts_DE[features_biotype, ep_samples, drop = FALSE]

  ## Remove features absent from all epididymis samples
  counts_bio_ep <- counts_bio_ep[
    rowSums(counts_bio_ep, na.rm = TRUE) > 0,
    ,
    drop = FALSE
  ]

  message("Features gardées après filtre >0 : ", nrow(counts_bio_ep))

  if (nrow(counts_bio_ep) < 2) {
    warning(
      "Biotype ", biotype_name,
      " ignoré après filtre : moins de 2 features exprimées dans l'épididyme."
    )
    return(NULL)
  }

  ## Avoid a DESeq2 error if a sample contains no counts
  samples_all_zero <- colSums(counts_bio_ep, na.rm = TRUE) == 0

  if (any(samples_all_zero)) {
    warning(
      "Biotype ", biotype_name,
      " ignoré : aucun count dans les échantillons ",
      paste(colnames(counts_bio_ep)[samples_all_zero], collapse = ", "),
      "."
    )
    return(NULL)
  }

  ## VST transformation
  dds_bio_ep <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts_bio_ep),
    colData = meta_ep_pca,
    design = ~ 1
  )

  vsd_bio_ep <- safe_vst(
    dds_bio_ep,
    label = paste0("PCA épididyme — biotype ", biotype_name)
  )

  vst_bio_ep <- SummarizedExperiment::assay(vsd_bio_ep)

  ## PCA
  pca_bio_ep <- prcomp(t(vst_bio_ep), scale. = FALSE)
  percentVar <- round(
    100 * (pca_bio_ep$sdev^2 / sum(pca_bio_ep$sdev^2)),
    1
  )

  pca_bio_ep_df <- as.data.frame(pca_bio_ep$x[, 1:2, drop = FALSE]) %>%
    tibble::rownames_to_column("A_sample") %>%
    dplyr::left_join(
      meta_ep_pca %>% as.data.frame(),
      by = "A_sample"
    )

  ## Save coordinates
  write.table(
    pca_bio_ep_df,
    file = file.path(
      pca_biotype_ep_dir,
      paste0("PCA_coordinates_Epididymis_", biotype_name, ".tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  ## Figure
  pdf(
    file.path(
      pca_biotype_ep_dir,
      paste0("PCA_Epididymis_", biotype_name, ".pdf")
    ),
    width = 9,
    height = 7,
    useDingbats = FALSE
  )

  p <- ggplot(
    pca_bio_ep_df,
    aes(
      x = PC1,
      y = PC2,
      color = `Origin in Epididymis`,
      label = `Unique Sample name`
    )
  ) +
    geom_point(size = 4, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 100) +
    scale_color_manual(values = cols_region, drop = FALSE) +
    labs(
      title = paste0("PCA epididymis — biotype ", biotype_name),
      subtitle = paste0(nrow(counts_bio_ep), " features utilisées"),
      x = paste0("PC1 : ", percentVar[1], "%"),
      y = paste0("PC2 : ", percentVar[2], "%"),
      color = "Origin in Epididymis"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )

  print(p)
  invisible(dev.off())

  return(pca_bio_ep_df)
}

## Keep only biotypes with at least one expressed feature
## in epididymis samples
features_present_ep <- rownames(counts_DE)[
  rowSums(counts_DE[, ep_samples, drop = FALSE], na.rm = TRUE) > 0
]

biotypes_to_pca_ep <- counts_filtered_out %>%
  dplyr::filter(
    feature_id %in% features_present_ep,
    !is.na(biotype_priority),
    biotype_priority != ""
  ) %>%
  dplyr::pull(biotype_priority) %>%
  unique()

print(biotypes_to_pca_ep)

pca_by_biotype_ep_results <- lapply(
  biotypes_to_pca_ep,
  make_pca_one_biotype_epididymis
)

names(pca_by_biotype_ep_results) <- biotypes_to_pca_ep

message("PCA par biotype dans l'épididyme terminée.")


## ----figures-13-pca-by-biotype-in-epididymis, echo=FALSE, message=FALSE, warning=FALSE, results="asis"-------------------------------
fig_dir <- file.path(outdir, "PCA_by_biotype_Epididymis")
figs <- if (dir.exists(fig_dir)) {
  list.files(fig_dir, pattern = "\\.pdf$", full.names = TRUE)
} else {
  character(0)
}
include_pdf_as_large_png(figs)


## ----14-pca-by-biotype-spz-and-epididymis, eval=TRUE---------------------------------------------------------------------------------
############################################################
# 14. PCA by biotype combining spermatozoa and epididymis
############################################################

message("---- PCA par biotype : spermatozoïdes + épididyme ----")

pca_biotype_spz_ep_dir <- file.path(
  outdir,
  "PCA_by_biotype_SPZ_Epididymis"
)
dir.create(
  pca_biotype_spz_ep_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

## Spermatozoa + epididymis metadata: Head / Body / Tail
meta_spz_ep_pca <- meta %>%
  dplyr::filter(
    `Tissue type` %in% c("Spermatozoa", "Epididymis"),
    `Origin in Epididymis` %in% c("Head", "Body", "Tail")
  )

meta_spz_ep_pca <- as.data.frame(meta_spz_ep_pca)
rownames(meta_spz_ep_pca) <- meta_spz_ep_pca$A_sample

## Samples present in the read-level matrix
spz_ep_samples <- intersect(
  rownames(meta_spz_ep_pca),
  colnames(counts_DE)
)
meta_spz_ep_pca <- meta_spz_ep_pca[spz_ep_samples, , drop = FALSE]

## Explicit ordering for legends
meta_spz_ep_pca$`Origin in Epididymis` <- factor(
  meta_spz_ep_pca$`Origin in Epididymis`,
  levels = c("Head", "Body", "Tail")
)

meta_spz_ep_pca$`Tissue type` <- factor(
  meta_spz_ep_pca$`Tissue type`,
  levels = c("Spermatozoa", "Epididymis")
)

message(
  "Nombre de samples SPZ + épididyme pour PCA : ",
  length(spz_ep_samples)
)

message("Répartition des samples par tissu et région :")
print(
  table(
    meta_spz_ep_pca$`Tissue type`,
    meta_spz_ep_pca$`Origin in Epididymis`,
    useNA = "ifany"
  )
)

if (length(spz_ep_samples) < 3) {
  stop("Moins de 3 échantillons SPZ + épididyme disponibles : PCA impossible.")
}

## Shared PCA function for SPZ + epididymis for one biotype
make_pca_one_biotype_spz_ep <- function(biotype_name) {

  message("PCA SPZ + épididyme — biotype : ", biotype_name)

  ## Features belonging to the requested biotype
  features_biotype <- counts_filtered_out %>%
    dplyr::filter(biotype_priority == biotype_name) %>%
    dplyr::pull(feature_id) %>%
    intersect(rownames(counts_DE))

  message(
    "Nombre de features pour ",
    biotype_name,
    " : ",
    length(features_biotype)
  )

  if (length(features_biotype) < 2) {
    warning("Biotype ", biotype_name, " ignoré : moins de 2 features.")
    return(NULL)
  }

  ## Biotype matrix restricted to SPZ + epididymis samples
  counts_bio_spz_ep <- counts_DE[
    features_biotype,
    spz_ep_samples,
    drop = FALSE
  ]

  ## Remove features absent from all retained samples
  counts_bio_spz_ep <- counts_bio_spz_ep[
    rowSums(counts_bio_spz_ep, na.rm = TRUE) > 0,
    ,
    drop = FALSE
  ]

  message(
    "Features gardées après filtre >0 : ",
    nrow(counts_bio_spz_ep)
  )

  if (nrow(counts_bio_spz_ep) < 2) {
    warning(
      "Biotype ", biotype_name,
      " ignoré après filtre : moins de 2 features exprimées."
    )
    return(NULL)
  }

  ## Avoid a normalization error if a sample contains no counts
  samples_all_zero <- colSums(counts_bio_spz_ep, na.rm = TRUE) == 0

  if (any(samples_all_zero)) {
    warning(
      "Biotype ", biotype_name,
      " ignoré : aucun count dans les échantillons ",
      paste(colnames(counts_bio_spz_ep)[samples_all_zero], collapse = ", "),
      "."
    )
    return(NULL)
  }

  ## VST transformation
  dds_bio_spz_ep <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts_bio_spz_ep),
    colData = meta_spz_ep_pca,
    design = ~ 1
  )

  vsd_bio_spz_ep <- safe_vst(
    dds_bio_spz_ep,
    label = paste0("PCA SPZ + épididyme — biotype ", biotype_name)
  )

  vst_bio_spz_ep <- SummarizedExperiment::assay(vsd_bio_spz_ep)

  ## Joint PCA of both tissues
  pca_bio_spz_ep <- prcomp(t(vst_bio_spz_ep), scale. = FALSE)
  percentVar <- round(
    100 * (pca_bio_spz_ep$sdev^2 / sum(pca_bio_spz_ep$sdev^2)),
    1
  )

  pca_bio_spz_ep_df <- as.data.frame(
    pca_bio_spz_ep$x[, 1:2, drop = FALSE]
  ) %>%
    tibble::rownames_to_column("A_sample") %>%
    dplyr::left_join(
      meta_spz_ep_pca %>% as.data.frame(),
      by = "A_sample"
    )

  ## Save coordinates
  write.table(
    pca_bio_spz_ep_df,
    file = file.path(
      pca_biotype_spz_ep_dir,
      paste0(
        "PCA_coordinates_SPZ_Epididymis_",
        biotype_name,
        ".tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  ## Figure
  pdf(
    file.path(
      pca_biotype_spz_ep_dir,
      paste0("PCA_SPZ_Epididymis_", biotype_name, ".pdf")
    ),
    width = 9,
    height = 7,
    useDingbats = FALSE
  )

  p <- ggplot(
    pca_bio_spz_ep_df,
    aes(
      x = PC1,
      y = PC2,
      color = `Origin in Epididymis`,
      shape = `Tissue type`,
      label = `Unique Sample name`
    )
  ) +
    geom_point(size = 4, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 100) +
    scale_color_manual(values = cols_region, drop = FALSE) +
    scale_shape_manual(
      values = c(Spermatozoa = 16, Epididymis = 17),
      drop = FALSE
    ) +
    labs(
      title = paste0(
        "PCA SPZ + epididymis — biotype ",
        biotype_name
      ),
      subtitle = paste0(
        nrow(counts_bio_spz_ep),
        " features utilisées"
      ),
      x = paste0("PC1 : ", percentVar[1], "%"),
      y = paste0("PC2 : ", percentVar[2], "%"),
      color = "Origin in Epididymis",
      shape = "Tissue type"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )

  print(p)
  invisible(dev.off())

  return(pca_bio_spz_ep_df)
}

## Keep biotypes with at least one expressed feature
## in SPZ + epididymis samples
features_present_spz_ep <- rownames(counts_DE)[
  rowSums(
    counts_DE[, spz_ep_samples, drop = FALSE],
    na.rm = TRUE
  ) > 0
]

biotypes_to_pca_spz_ep <- counts_filtered_out %>%
  dplyr::filter(
    feature_id %in% features_present_spz_ep,
    !is.na(biotype_priority),
    biotype_priority != ""
  ) %>%
  dplyr::pull(biotype_priority) %>%
  unique()

print(biotypes_to_pca_spz_ep)

pca_by_biotype_spz_ep_results <- lapply(
  biotypes_to_pca_spz_ep,
  make_pca_one_biotype_spz_ep
)

names(pca_by_biotype_spz_ep_results) <- biotypes_to_pca_spz_ep

message("PCA par biotype SPZ + épididyme terminée.")


## ----figures-14-pca-by-biotype-spz-and-epididymis, echo=FALSE, message=FALSE, warning=FALSE, results="asis"--------------------------
fig_dir <- file.path(outdir, "PCA_by_biotype_SPZ_Epididymis")
figs <- if (dir.exists(fig_dir)) {
  list.files(fig_dir, pattern = "\\.pdf$", full.names = TRUE)
} else {
  character(0)
}
include_pdf_as_large_png(figs)


## ----15-read-level-deseq2-spz-head-vs-tail-cdna-volcano, eval=TRUE-------------------------------------------------------------------
############################################################
# 15. Read-level DESeq2 SPZ Head vs Tail + cDNA volcano
############################################################

message("---- Read-level DESeq2 SPZ Head vs Tail ----")

read_volcano_dir <- file.path(outdir, "Volcano_read_level_SPZ")
dir.create(read_volcano_dir, showWarnings = FALSE, recursive = TRUE)

meta_spz_ht <- prepare_spz_metadata(regions = c("Head", "Tail"))
meta_spz_ht$Origin_in_Epididymis <- relevel(meta_spz_ht$Origin_in_Epididymis, ref = "Head")

counts_spz_ht <- counts_DE[, rownames(meta_spz_ht), drop = FALSE]
stopifnot(all(colnames(counts_spz_ht) == rownames(meta_spz_ht)))

dds_spz_ht <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(counts_spz_ht),
  colData = meta_spz_ht,
  design = ~ Origin_in_Epididymis
)

dds_spz_ht <- DESeq2::DESeq(dds_spz_ht)

res_head_vs_tail <- DESeq2::results(
  dds_spz_ht,
  contrast = c("Origin_in_Epididymis", "Head", "Tail")
)

res_head_vs_tail_df <- as.data.frame(res_head_vs_tail) %>%
  tibble::rownames_to_column("feature_id")

write.table(
  res_head_vs_tail_df,
  file = file.path(read_volcano_dir, "DESeq2_results_read_level_SPZ_Head_vs_Tail.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

annot_reads <- counts_filtered_out %>%
  dplyr::select(feature_id, sequence, gene_name_priority, biotype_priority)

res_head_vs_tail_annot <- res_head_vs_tail_df %>%
  dplyr::left_join(annot_reads, by = "feature_id") %>%
  dplyr::mutate(
    contrast = "Head_vs_Tail_read_level",
    neglog10padj = -log10(padj),

    ## The DESeq2 contrast remains Head vs Tail:
    ## log2FoldChange > 0 = Head > Tail
    ## log2FoldChange < 0 = Tail > Head
    ## Only the sign is inverted for volcano plot display so that
    ## Tail appears on the right, without modifying the DESeq2 results.
    plot_log2FoldChange = -log2FoldChange,

    signif = !is.na(padj) & padj < padj_cut & abs(log2FoldChange) >= lfc_cut,
    direction = dplyr::case_when(
      signif & log2FoldChange > 0 ~ "Higher_in_Head",
      signif & log2FoldChange < 0 ~ "Higher_in_Tail",
      TRUE ~ "Not_DE"
    )
  )

write.table(
  res_head_vs_tail_annot,
  file = file.path(read_volcano_dir, "Volcano_table_read_level_SPZ_Head_vs_Tail.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## Read-level volcano plot for all biotypes
volc_head_vs_tail <- res_head_vs_tail_annot

n_de_read <- sum(volc_head_vs_tail$signif, na.rm = TRUE)

## Total number of features included in the volcano plot
n_input_read <- nrow(volc_head_vs_tail)

## Number of features with an available adjusted p-value
n_with_padj_read <- sum(!is.na(volc_head_vs_tail$padj))

top_labels_read <- volc_head_vs_tail %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 25)

pdf(file.path(read_volcano_dir, "Volcano_read_level_SPZ_Head_vs_Tail.pdf"),
    width = 10, height = 9, useDingbats = FALSE)

p_read <- ggplot(
  volc_head_vs_tail,
  aes(
    x = plot_log2FoldChange,
    y = neglog10padj,
    color = direction
  )
) +
  geom_point(
    alpha = 0.65,
    size = 1.5
  ) +
  scale_color_manual(
    values = c(
      "Not_DE" = "#C7C7C7",
      "Higher_in_Head" = "#617C93",
      "Higher_in_Tail" = "#A77768"
    )
  ) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.6) +
  ggrepel::geom_text_repel(
    data = top_labels_read,
    aes(label = feature_id),
    size = 3.0,
    max.overlaps = 50,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0
  ) +
  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = paste0(
      "DE features: ", n_de_read,
      "\nInput features: ", n_input_read,
      "\nFeatures with padj: ", n_with_padj_read
    ),
    hjust = 1.02,
    vjust = 1.02,
    size = 3.6,
    label.size = 0.25,
    fill = "white",
    alpha = 0.9,
    lineheight = 0.9
  ) +
  labs(
    title = "Volcano read-level SPZ: Head vs Tail",
    subtitle = paste0("feature_id = sequence + gene | padj < ", padj_cut, " et |log2FC| ≥ ", lfc_cut),
    x = "Displayed log2(FoldChange): Tail / Head",
    y = "-log10(padj)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())

print(p_read)
invisible(dev.off())

## cDNA read-level only
volc_cdna_head_vs_tail <- res_head_vs_tail_annot %>%
  dplyr::filter(biotype_priority == "cdna")

n_cdna_de <- sum(volc_cdna_head_vs_tail$signif, na.rm = TRUE)

## Total number of cDNA features included in the volcano plot
n_input_cdna <- nrow(volc_cdna_head_vs_tail)

## Number of cDNA features with an available adjusted p-value
n_with_padj_cdna <- sum(!is.na(volc_cdna_head_vs_tail$padj))

message(
  "Nombre de cDNA features DE : ", n_cdna_de,
  " | Input cDNA features : ", n_input_cdna,
  " | cDNA features avec padj : ", n_with_padj_cdna
)

cdna_de_direction_summary <- volc_cdna_head_vs_tail %>%
  dplyr::filter(signif) %>%
  dplyr::count(direction, name = "n_cdna_DE_features")

write.table(
  volc_cdna_head_vs_tail,
  file = file.path(read_volcano_dir, "Volcano_table_cdna_reads_SPZ_Head_vs_Tail_read_level.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  cdna_de_direction_summary,
  file = file.path(read_volcano_dir, "Volcano_cdna_reads_SPZ_Head_vs_Tail_direction_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

top_cdna_labels <- volc_cdna_head_vs_tail %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = topN)

pdf(file.path(read_volcano_dir, "Volcano_plot_cdna_reads_SPZ_Head_vs_Tail_read_level.pdf"),
    width = 10, height = 9, useDingbats = FALSE)

p_cdna <- ggplot(
  volc_cdna_head_vs_tail,
  aes(
    x = plot_log2FoldChange,
    y = neglog10padj,
    color = direction
  )
) +
  geom_point(
    alpha = 0.65,
    size = 1.5
  ) +
  scale_color_manual(
    values = c(
      "Not_DE" = "#C7C7C7",
      "Higher_in_Head" = "#617C93",
      "Higher_in_Tail" = "#A77768"
    )
  ) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.6) +
  ggrepel::geom_text_repel(
    data = top_cdna_labels,
    aes(label = feature_id),
    size = 3.0,
    max.overlaps = 50,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0
  ) +
  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = paste0(
      "DE cDNA features: ", n_cdna_de,
      "\nInput cDNA features: ", n_input_cdna,
      "\ncDNA features with padj: ", n_with_padj_cdna
    ),
    hjust = 1.02,
    vjust = 1.02,
    size = 3.6,
    label.size = 0.25,
    fill = "white",
    alpha = 0.9,
    lineheight = 0.9
  ) +
  labs(
    title = "Volcano: cDNA reads/features — Head vs Tail (SPZ)",
    subtitle = paste0("Read-level | padj < ", padj_cut, " et |log2FC| ≥ ", lfc_cut),
    x = "Displayed log2(FoldChange): Tail / Head",
    y = "-log10(padj)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())

print(p_cdna)
invisible(dev.off())

## Diagnostic of the numbers displayed in the read-level volcano plots
read_level_head_vs_tail_summary <- tibble::tibble(
  contrast = "Head_vs_Tail",
  level = c("all_read_features", "cdna_read_features"),
  n_input_features = c(n_input_read, n_input_cdna),
  n_with_pvalue = c(
    sum(!is.na(volc_head_vs_tail$pvalue)),
    sum(!is.na(volc_cdna_head_vs_tail$pvalue))
  ),
  n_with_padj = c(n_with_padj_read, n_with_padj_cdna),
  n_DE_features = c(n_de_read, n_cdna_de)
)

print(read_level_head_vs_tail_summary)

write.table(
  read_level_head_vs_tail_summary,
  file = file.path(read_volcano_dir, "diagnostic_read_level_Head_vs_Tail_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## ----figures-15-read-level-deseq2-spz-head-vs-tail-cdna-volcano, echo=FALSE, message=FALSE, warning=FALSE, results="asis"------------
figs <- c(
  file.path(outdir, "Volcano_read_level_SPZ/Volcano_read_level_SPZ_Head_vs_Tail.pdf"),
  file.path(outdir, "Volcano_read_level_SPZ/Volcano_plot_cdna_reads_SPZ_Head_vs_Tail_read_level.pdf")
)
include_pdf_as_large_png(figs)


## ----16-gene-level-deseq2-spz-head-vs-body-head-vs-tail-body-vs-t, eval=TRUE---------------------------------------------------------
############################################################
# 16. Gene-level DESeq2 SPZ: Head vs Body, Head vs Tail, Body vs Tail
############################################################

message("---- Gene-level DESeq2 SPZ toutes comparaisons ----")

gene_volcano_dir <- file.path(outdir, "Volcano_gene_level_SPZ_all_genes")
dir.create(gene_volcano_dir, showWarnings = FALSE, recursive = TRUE)

meta_spz_all_regions <- prepare_spz_metadata(regions = c("Head", "Body", "Tail"))

message("Samples SPZ par région :")
print(table(meta_spz_all_regions$Origin_in_Epididymis, useNA = "ifany"))

counts_gene_spz_all_regions <- counts_gene_DE[, rownames(meta_spz_all_regions), drop = FALSE]
stopifnot(all(colnames(counts_gene_spz_all_regions) == rownames(meta_spz_all_regions)))

dds_gene_spz_regions <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(counts_gene_spz_all_regions),
  colData = meta_spz_all_regions,
  design = ~ Origin_in_Epididymis
)

dds_gene_spz_regions <- DESeq2::DESeq(dds_gene_spz_regions)

get_gene_contrast <- function(group1, group2) {

  contrast_name <- paste0(group1, "_vs_", group2)
  message("Contraste gene-level : ", contrast_name)

  res <- DESeq2::results(
    dds_gene_spz_regions,
    contrast = c("Origin_in_Epididymis", group1, group2)
  )

  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(annot_gene, by = "gene_id") %>%
    dplyr::mutate(contrast = contrast_name)

  write.table(
    res_df,
    file = file.path(gene_volcano_dir, paste0("DESeq2_gene_level_SPZ_", contrast_name, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  return(res_df)
}



plot_gene_volcano <- function(res_df, contrast_name, group1, group2, out_dir, biotype_name = NULL) {
  volcano_colors <- c(
    "Not_DE" = "#C7C7C7",
    "#617C93",
    "#A77768"
  )
  
  names(volcano_colors)[2:3] <- c(
    paste0("Higher_in_", group1),
    paste0("Higher_in_", group2)
  )
  volc <- res_df %>%
    dplyr::mutate(
      neglog10padj = -log10(padj),

      ## The DESeq2 contrast remains group1 vs group2.
      ## The original log2FoldChange sign is retained in the table.
      ## This variable is used only to display the second group
      ## on the right side of the volcano plot (Body for HB, Tail for BT and HT).
      plot_log2FoldChange = -log2FoldChange,
      
      signif = !is.na(padj) & padj < padj_cut & abs(log2FoldChange) >= lfc_cut,
      direction = dplyr::case_when(
        signif & log2FoldChange > 0 ~ paste0("Higher_in_", group1),
        signif & log2FoldChange < 0 ~ paste0("Higher_in_", group2),
        TRUE ~ "Not_DE"
      )
    )

  if (!is.null(biotype_name)) {
    volc <- volc %>% dplyr::filter(biotype_priority == biotype_name)
  }

  if (nrow(volc) == 0) {
    warning("Aucune ligne pour ", contrast_name, " ", biotype_name)
    return(NULL)
  }

  n_de <- sum(volc$signif, na.rm = TRUE)
  n_input <- nrow(volc)
  n_with_padj <- sum(!is.na(volc$padj))

  label_suffix <- ifelse(is.null(biotype_name), "all_genes", biotype_name)

  message(
    contrast_name, " — ", label_suffix,
    " : DE genes ", n_de,
    " | input genes ", n_input,
    " | genes with padj ", n_with_padj
  )

  write.table(
    volc,
    file = file.path(out_dir, paste0("Volcano_table_gene_level_SPZ_", contrast_name, "_", label_suffix, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  de_summary <- volc %>%
    dplyr::filter(signif) %>%
    dplyr::count(direction, name = "n_DE_genes")

  write.table(
    de_summary,
    file = file.path(out_dir, paste0("DE_summary_gene_level_SPZ_", contrast_name, "_", label_suffix, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  top_labels <- volc %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = topN)

  title_extra <- ifelse(is.null(biotype_name), "all genes", biotype_name)

  pdf(file.path(out_dir, paste0("Volcano_gene_level_SPZ_", contrast_name, "_", label_suffix, ".pdf")),
      width = 10, height = 9, useDingbats = FALSE)

  p <- ggplot(
    volc,
    aes(
      x = plot_log2FoldChange,
      y = neglog10padj,
      color = direction
    )
  ) +
    geom_point(
      alpha = 0.65,
      size = 1.5
    ) +
    scale_color_manual(
      values = volcano_colors,
      breaks = c(
        "Not_DE",
        paste0("Higher_in_", group1),
        paste0("Higher_in_", group2)
      )
    ) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.6) +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.6) +
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(label = gene_name_priority),
      size = 3.2,
      max.overlaps = 50,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0
    ) +
    annotate(
      "label",
      x = Inf,
      y = Inf,
      label = paste0(
        "DE genes: ", n_de,
        "\nInput genes: ", n_input,
        "\nGenes with padj: ", n_with_padj
      ),
      hjust = 1.02,
      vjust = 1.02,
      size = 3.6,
      label.size = 0.25,
      fill = "white",
      alpha = 0.9,
      lineheight = 0.9
    ) +
    labs(
      title = paste0("Volcano gene-level SPZ: ", group1, " vs ", group2, " — ", title_extra),
      subtitle = paste0(
        "DESeq2 contrast = ", group1, " vs ", group2,
        " | x-axis inverted for display | padj < ", padj_cut,
        " et |log2FC| ≥ ", lfc_cut
      ),
      x = paste0("Displayed log2(FoldChange): ", group2, " / ", group1),
      y = "-log10(padj)"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())

  print(p)
  invisible(dev.off())

  return(volc)
}

comparisons_spz <- list(
  c("Head", "Body"),
  c("Head", "Tail"),
  c("Body", "Tail")
)

gene_level_volcano_results <- list()

for (comp in comparisons_spz) {

  group1 <- comp[1]
  group2 <- comp[2]
  contrast_name <- paste0(group1, "_vs_", group2)

  res_df <- get_gene_contrast(group1, group2)
  volc_df <- plot_gene_volcano(
    res_df = res_df,
    contrast_name = contrast_name,
    group1 = group1,
    group2 = group2,
    out_dir = gene_volcano_dir,
    biotype_name = NULL
  )

  gene_level_volcano_results[[contrast_name]] <- volc_df
}

message("Volcano plots gene-level toutes comparaisons terminés.")

## Diagnostic: the total number of input genes must be identical
## across gene-level all-gene contrasts. The n_with_padj column may vary
## because DESeq2 may assign padj = NA to some genes depending on the contrast.
gene_level_contrast_summary <- lapply(
  names(gene_level_volcano_results),
  function(contrast_name) {
    df <- gene_level_volcano_results[[contrast_name]]
    tibble::tibble(
      contrast = contrast_name,
      n_input_genes = nrow(df),
      n_with_pvalue = sum(!is.na(df$pvalue)),
      n_with_padj = sum(!is.na(df$padj)),
      n_DE_genes = sum(
        !is.na(df$padj) &
          df$padj < padj_cut &
          abs(df$log2FoldChange) >= lfc_cut
      )
    )
  }
) %>%
  dplyr::bind_rows()

print(gene_level_contrast_summary)

write.table(
  gene_level_contrast_summary,
  file = file.path(gene_volcano_dir, "diagnostic_gene_level_contrast_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## ----figures-16-gene-level-deseq2-spz-head-vs-body-head-vs-tail-body-vs-t, echo=FALSE, message=FALSE, warning=FALSE, results="asis"----
fig_dir <- file.path(outdir, "Volcano_gene_level_SPZ_all_genes")
figs <- if (dir.exists(fig_dir)) list.files(fig_dir, pattern = "\\.pdf$", full.names = TRUE) else character(0)
include_pdf_as_large_png(figs)


## ----17-gene-level-volcano-plots-by-biotype, eval=TRUE-------------------------------------------------------------------------------
############################################################
# 17. Gene-level volcano plots by biotype
############################################################

message("---- Volcano plots par biotype gene-level SPZ ----")

biotype_volcano_dir <- file.path(outdir, "Volcano_gene_level_SPZ_by_biotype")
dir.create(biotype_volcano_dir, showWarnings = FALSE, recursive = TRUE)

biotypes_for_volcano <- annot_gene %>%
  dplyr::filter(!is.na(biotype_priority), biotype_priority != "") %>%
  dplyr::pull(biotype_priority) %>%
  unique()

print(biotypes_for_volcano)

biotype_volcano_results <- list()

for (comp in comparisons_spz) {

  group1 <- comp[1]
  group2 <- comp[2]
  contrast_name <- paste0(group1, "_vs_", group2)
  res_df <- get_gene_contrast(group1, group2)

  for (biotype_name in biotypes_for_volcano) {

    key_name <- paste0(contrast_name, "_", biotype_name)

    biotype_volcano_results[[key_name]] <- plot_gene_volcano(
      res_df = res_df,
      contrast_name = contrast_name,
      group1 = group1,
      group2 = group2,
      out_dir = biotype_volcano_dir,
      biotype_name = biotype_name
    )
  }
}

de_summary_all_biotypes <- lapply(
  names(biotype_volcano_results),
  function(key_name) {

    df <- biotype_volcano_results[[key_name]]

    if (is.null(df)) {
      return(NULL)
    }

    df %>%
      dplyr::filter(signif) %>%
      dplyr::group_by(contrast, biotype_priority, direction) %>%
      dplyr::summarise(n_DE_genes = dplyr::n(), .groups = "drop")
  }
) %>%
  dplyr::bind_rows()

print(de_summary_all_biotypes)

write.table(
  de_summary_all_biotypes,
  file = file.path(biotype_volcano_dir, "DE_summary_all_comparisons_by_biotype.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Script terminé. Résultats dans : ", outdir)


## ----figures-17-gene-level-volcano-plots-by-biotype, echo=FALSE, message=FALSE, warning=FALSE, results="asis"------------------------
fig_dir <- file.path(outdir, "Volcano_gene_level_SPZ_by_biotype")
figs <- if (dir.exists(fig_dir)) list.files(fig_dir, pattern = "\\.pdf$", full.names = TRUE) else character(0)
include_pdf_as_large_png(figs)


## ----18-barplot-total-genes-and-up-down-de-by-biotype, eval=TRUE---------------------------------------------------------------------
############################################################
# 18. Barplot of total genes and up/downregulated DE genes
#     by biotype
############################################################

message(
  paste0(
    "---- Barplot du nombre total de gènes et des gènes DE ",
    "up/down par biotype ----"
  )
)

barplot_de_biotype_dir <- file.path(
  outdir,
  "Barplot_DE_genes_by_biotype"
)

dir.create(
  barplot_de_biotype_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

## Requested order for the three comparisons
contrast_order_barplot <- c(
  "Head_vs_Body",
  "Body_vs_Tail",
  "Head_vs_Tail"
)

comparison_labels_barplot <- c(
  Head_vs_Body = "HB",
  Body_vs_Tail = "BT",
  Head_vs_Tail = "HT"
)

## Logical biotype order according to annotation priority
preferred_biotype_order <- c(
  "miR",
  "tRF",
  "rRF",
  "snoRNA",
  "piR",
  "ncRNA",
  "cdna"
)

observed_biotypes_gene <- annot_gene %>%
  dplyr::filter(
    !is.na(biotype_priority),
    biotype_priority != "",
    !is.na(gene_id),
    gene_id != ""
  ) %>%
  dplyr::pull(biotype_priority) %>%
  unique() %>%
  as.character()

biotype_order_barplot <- c(
  intersect(preferred_biotype_order, observed_biotypes_gene),
  setdiff(sort(observed_biotypes_gene), preferred_biotype_order)
)

## Total number of genes belonging to each biotype.
## This number does not depend on the contrast.
total_genes_by_biotype <- annot_gene %>%
  dplyr::filter(
    !is.na(biotype_priority),
    biotype_priority != "",
    !is.na(gene_id),
    gene_id != ""
  ) %>%
  dplyr::distinct(gene_id, biotype_priority) %>%
  dplyr::count(
    biotype_priority,
    name = "n_total_genes"
  ) %>%
  dplyr::mutate(
    biotype_priority = as.character(biotype_priority)
  )

print(total_genes_by_biotype)

## Number of up- and downregulated DE genes for each contrast.
## The signif object already includes both thresholds:
## padj < padj_cut and abs(log2FoldChange) >= lfc_cut.
## DESeq2 contrasts remain unchanged (group1 vs group2).
## To harmonize the figures with the Head → Body → Tail progression:
##   UP   = higher expression in the second group = log2FoldChange < 0
##   DOWN = higher expression in the first group = log2FoldChange > 0
## No DESeq2 result or matrix is modified.
de_counts_by_contrast <- lapply(
  contrast_order_barplot,
  function(contrast_name) {

    df <- gene_level_volcano_results[[contrast_name]]

    if (is.null(df)) {
      stop(
        "Résultats gene-level absents pour le contraste : ",
        contrast_name
      )
    }

    df %>%
      dplyr::filter(
        !is.na(biotype_priority),
        biotype_priority != "",
        !is.na(gene_id),
        gene_id != ""
      ) %>%
      dplyr::group_by(biotype_priority) %>%
      dplyr::summarise(
        n_DE_total = dplyr::n_distinct(
          gene_id[signif %in% TRUE]
        ),
        n_up_DE = dplyr::n_distinct(
          gene_id[
            signif %in% TRUE &
              !is.na(log2FoldChange) &
              log2FoldChange < 0
          ]
        ),
        n_down_DE = dplyr::n_distinct(
          gene_id[
            signif %in% TRUE &
              !is.na(log2FoldChange) &
              log2FoldChange > 0
          ]
        ),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        contrast = contrast_name,
        biotype_priority = as.character(biotype_priority)
      )
  }
) %>%
  dplyr::bind_rows()

## Build all biotype x contrast combinations,
## including those containing no DE genes.
de_barplot_summary <- tidyr::expand_grid(
  contrast = contrast_order_barplot,
  biotype_priority = biotype_order_barplot
) %>%
  dplyr::left_join(
    total_genes_by_biotype,
    by = "biotype_priority"
  ) %>%
  dplyr::left_join(
    de_counts_by_contrast,
    by = c("contrast", "biotype_priority")
  ) %>%
  dplyr::mutate(
    n_total_genes = tidyr::replace_na(n_total_genes, 0L),
    n_DE_total = tidyr::replace_na(n_DE_total, 0L),
    n_up_DE = tidyr::replace_na(n_up_DE, 0L),
    n_down_DE = tidyr::replace_na(n_down_DE, 0L),
    comparison = factor(
      unname(comparison_labels_barplot[contrast]),
      levels = c("HB", "BT", "HT")
    ),
    biotype_priority = factor(
      biotype_priority,
      levels = biotype_order_barplot
    )
  ) %>%
  dplyr::arrange(biotype_priority, comparison)

## Consistency checks before plotting.
if (any(
  de_barplot_summary$n_DE_total !=
    de_barplot_summary$n_up_DE + de_barplot_summary$n_down_DE
)) {
  stop(
    "Incohérence dans le bloc 18 : n_DE_total doit être égal ",
    "à n_up_DE + n_down_DE."
  )
}

if (any(
  de_barplot_summary$n_DE_total >
    de_barplot_summary$n_total_genes
)) {
  stop(
    "Incohérence dans le bloc 18 : le nombre de gènes DE ",
    "dépasse le nombre total de gènes d'un biotype."
  )
}

## Other genes correspond to non-DE genes under the thresholds used.
de_barplot_summary <- de_barplot_summary %>%
  dplyr::mutate(
    n_other_genes = n_total_genes - n_DE_total
  )

print(de_barplot_summary)

write.table(
  de_barplot_summary,
  file = file.path(
    barplot_de_biotype_dir,
    paste0(
      "Total_genes_and_up_down_DE_genes_",
      "by_biotype_and_comparison.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## Convert to long format to build a single stacked bar:
## other genes + downregulated DE genes + upregulated DE genes = total genes.
de_barplot_long <- de_barplot_summary %>%
  dplyr::select(
    contrast,
    comparison,
    biotype_priority,
    n_total_genes,
    n_up_DE,
    n_down_DE,
    n_other_genes
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      n_other_genes,
      n_down_DE,
      n_up_DE
    ),
    names_to = "gene_status",
    values_to = "n_genes"
  ) %>%
  dplyr::mutate(
    gene_status = factor(
      gene_status,
      levels = c(
        "n_other_genes",
        "n_down_DE",
        "n_up_DE"
      ),
      labels = c(
        "Other genes",
        "Downregulated DE genes",
        "Upregulated DE genes"
      )
    )
  )

## More readable names for facet labels
biotype_facet_labels <- c(
  miR = "miRNA",
  tRF = "tRF",
  rRF = "rRF",
  snoRNA = "snoRNA",
  piR = "piRNA",
  ncRNA = "ncRNA",
  cdna = "cDNA"
)

pdf(
  file.path(
    barplot_de_biotype_dir,
    "Barplot_total_genes_and_up_down_DE_by_biotype.pdf"
  ),
  width = 13,
  height = 9,
  useDingbats = FALSE
)

p_de_barplot_biotype <- ggplot(
  de_barplot_long,
  aes(
    x = comparison,
    y = n_genes,
    fill = gene_status
  )
) +
  ## A single stacked bar whose height equals
  ## the total number of genes in the biotype.
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.35
  ) +
  ## Display the total and up/down counts above each bar.
  geom_text(
    data = de_barplot_summary %>%
      dplyr::filter(n_total_genes > 0),
    aes(
      x = comparison,
      y = n_total_genes,
      label = paste0(
        "Total: ", n_total_genes,
        "\nUp: ", n_up_DE,
        " | Down: ", n_down_DE
      )
    ),
    inherit.aes = FALSE,
    vjust = -0.25,
    size = 3.0,
    lineheight = 0.9
  ) +
  facet_wrap(
    ~ biotype_priority,
    scales = "free_y",
    ncol = 3,
    drop = FALSE,
    labeller = ggplot2::as_labeller(
      biotype_facet_labels,
      default = ggplot2::label_value
    )
  ) +
  scale_fill_manual(
    values = c(
      "Other genes" = "grey84",
      "Downregulated DE genes" = "#4C6A92",
      "Upregulated DE genes" = "#A66A5B"
    ),
    breaks = c(
      "Other genes",
      "Upregulated DE genes",
      "Downregulated DE genes"
    )
  ) +
  scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.28))
  ) +
  labs(
    title = paste0(
      "Total genes and up/downregulated DE genes ",
      "by biotype"
    ),
    subtitle = paste0(
      "HB = Head vs Body; BT = Body vs Tail; HT = Head vs Tail | ",
      "Up: higher in the second group; Down: higher in the first group"
    ),
    x = "Comparison",
    y = "Number of genes",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

print(p_de_barplot_biotype)
invisible(dev.off())

message(
  paste0(
    "Barplot du nombre total de gènes et des gènes DE ",
    "up/down terminé."
  )
)
message(
  "Figure générée : ",
  file.path(
    barplot_de_biotype_dir,
    "Barplot_total_genes_and_up_down_DE_by_biotype.pdf"
  )
)


## ----figures-18-barplot-total-genes-and-up-down-de-by-biotype, echo=FALSE, message=FALSE, warning=FALSE, results="asis"--------------
figs <- c(
  file.path(
    outdir,
    "Barplot_DE_genes_by_biotype",
    "Barplot_total_genes_and_up_down_DE_by_biotype.pdf"
  )
)
include_pdf_as_large_png(figs)


## ----19-venn-up-genes-by-biotype, eval=TRUE, message=FALSE, warning=FALSE------------------------------------------------------------
############################################################
# 19. Venn diagrams of upregulated DE genes by biotype
############################################################

message("---- Diagrammes de Venn des gènes DE up par biotype ----")

## The VennDiagram package is installed in block 1.
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  stop("Le package VennDiagram n'est pas disponible après le bloc d'installation.")
}

venn_up_dir <- file.path(
  outdir,
  "Venn_gene_level_SPZ_up_by_biotype"
)

dir.create(
  venn_up_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

## Remove old figures to avoid displaying obsolete files
old_venn_files <- list.files(
  venn_up_dir,
  pattern = "^Venn_up_genes_.*\\.pdf$",
  full.names = TRUE
)

if (length(old_venn_files) > 0) {
  invisible(file.remove(old_venn_files))
}

if (!exists("gene_level_volcano_results")) {
  stop(
    "L'objet gene_level_volcano_results est introuvable. ",
    "Le bloc gene-level DESeq2 doit être exécuté avant ce bloc."
  )
}

contrast_map_venn <- c(
  HB = "Head_vs_Body",
  BT = "Body_vs_Tail",
  HT = "Head_vs_Tail"
)

missing_contrasts_venn <- setdiff(
  unname(contrast_map_venn),
  names(gene_level_volcano_results)
)

if (length(missing_contrasts_venn) > 0) {
  stop(
    "Contrastes manquants dans gene_level_volcano_results : ",
    paste(missing_contrasts_venn, collapse = ", ")
  )
}

## Retrieve upregulated DE genes in the second group for a contrast and biotype
get_up_genes_for_venn <- function(contrast_name, biotype_name) {

  contrast_df <- gene_level_volcano_results[[contrast_name]]

  if (is.null(contrast_df) || nrow(contrast_df) == 0) {
    return(character(0))
  }

  contrast_df %>%
    dplyr::filter(
      as.character(biotype_priority) == biotype_name,
      signif %in% TRUE,
      !is.na(log2FoldChange),
      log2FoldChange < 0,
      !is.na(gene_id),
      gene_id != ""
    ) %>%
    dplyr::distinct(gene_id) %>%
    dplyr::pull(gene_id) %>%
    unique()
}

## Logical biotype order
preferred_biotype_order_venn <- c(
  "miR",
  "tRF",
  "rRF",
  "snoRNA",
  "piR",
  "ncRNA",
  "cdna"
)

observed_biotypes_venn <- annot_gene %>%
  dplyr::filter(
    !is.na(biotype_priority),
    biotype_priority != ""
  ) %>%
  dplyr::distinct(biotype_priority) %>%
  dplyr::pull(biotype_priority) %>%
  as.character()

biotypes_for_venn <- c(
  intersect(preferred_biotype_order_venn, observed_biotypes_venn),
  setdiff(sort(observed_biotypes_venn), preferred_biotype_order_venn)
)

print(biotypes_for_venn)

venn_up_summary_list <- list()
venn_up_membership_list <- list()

for (biotype_name in biotypes_for_venn) {

  message("Diagramme de Venn — biotype : ", biotype_name)

  venn_sets <- list(
    HB = get_up_genes_for_venn(
      contrast_name = contrast_map_venn[["HB"]],
      biotype_name = biotype_name
    ),
    BT = get_up_genes_for_venn(
      contrast_name = contrast_map_venn[["BT"]],
      biotype_name = biotype_name
    ),
    HT = get_up_genes_for_venn(
      contrast_name = contrast_map_venn[["HT"]],
      biotype_name = biotype_name
    )
  )

  ## Ensure unique character vectors without NA values
  venn_sets <- lapply(
    venn_sets,
    function(x) unique(as.character(x[!is.na(x) & x != ""]))
  )

  genes_HB <- venn_sets$HB
  genes_BT <- venn_sets$BT
  genes_HT <- venn_sets$HT

  genes_HB_BT <- intersect(genes_HB, genes_BT)
  genes_HB_HT <- intersect(genes_HB, genes_HT)
  genes_BT_HT <- intersect(genes_BT, genes_HT)
  genes_HB_BT_HT <- Reduce(
    intersect,
    list(genes_HB, genes_BT, genes_HT)
  )

  all_up_genes <- sort(
    unique(c(genes_HB, genes_BT, genes_HT))
  )

  venn_up_summary_list[[biotype_name]] <- tibble::tibble(
    biotype_priority = biotype_name,
    n_HB = length(genes_HB),
    n_BT = length(genes_BT),
    n_HT = length(genes_HT),
    n_HB_BT = length(genes_HB_BT),
    n_HB_HT = length(genes_HB_HT),
    n_BT_HT = length(genes_BT_HT),
    n_HB_BT_HT = length(genes_HB_BT_HT),
    n_union = length(all_up_genes)
  )

  if (length(all_up_genes) > 0) {

    membership_biotype <- tibble::tibble(
      gene_id = all_up_genes,
      HB = all_up_genes %in% genes_HB,
      BT = all_up_genes %in% genes_BT,
      HT = all_up_genes %in% genes_HT
    ) %>%
      dplyr::left_join(
        annot_gene %>%
          dplyr::select(
            gene_id,
            gene_name_priority,
            biotype_priority
          ) %>%
          dplyr::distinct(gene_id, .keep_all = TRUE),
        by = "gene_id"
      ) %>%
      dplyr::select(
        biotype_priority,
        gene_id,
        gene_name_priority,
        HB,
        BT,
        HT
      )

    venn_up_membership_list[[biotype_name]] <- membership_biotype
  }

  ## No diagram to generate if all three sets are empty
  if (length(all_up_genes) == 0) {
    message(
      "Aucun gène DE up pour le biotype ",
      biotype_name,
      " : aucun diagramme généré."
    )
    next
  }

  biotype_file_name <- gsub(
    "[^A-Za-z0-9_-]+",
    "_",
    biotype_name
  )

  venn_pdf_file <- file.path(
    venn_up_dir,
    paste0(
      "Venn_up_genes_",
      biotype_file_name,
      "_HB_BT_HT.pdf"
    )
  )

  ## filename = NULL returns the grid object without attempting to write directly
  ## a PDF format that is not supported by the imagetype argument.
  venn_grob <- VennDiagram::venn.diagram(
    x = venn_sets,
    filename = NULL,
    disable.logging = TRUE,
    category.names = c(
      "HB\nBody > Head",
      "BT\nTail > Body",
      "HT\nTail > Head"
    ),
    fill = c(
      "#F1D21A",
      "#F28E1C",
      "#D90D0D"
    ),
    alpha = c(0.45, 0.45, 0.45),
    col = c(
      "#8C7800",
      "#A95400",
      "#850000"
    ),
    lwd = c(2, 2, 2),
    cex = 1.3,
    fontface = "bold",
    cat.cex = c(1.1, 1.1, 1.1),
    cat.fontface = c("bold", "bold", "bold"),
    cat.col = c("black", "black", "black"),
    margin = 0.10,
    scaled = FALSE,
    euler.d = FALSE
  )

  grDevices::pdf(
    venn_pdf_file,
    width = 8,
    height = 8,
    useDingbats = FALSE
  )

  grid::grid.newpage()
  grid::grid.draw(venn_grob)
  grid::grid.text(
    paste0("Genes DE up toward the next region — biotype ", biotype_name),
    x = 0.5,
    y = 0.98,
    gp = grid::gpar(fontsize = 15, fontface = "bold")
  )

  invisible(grDevices::dev.off())
}

venn_up_summary <- dplyr::bind_rows(
  venn_up_summary_list
)

print(venn_up_summary)

write.table(
  venn_up_summary,
  file = file.path(
    venn_up_dir,
    "Venn_up_genes_summary_by_biotype.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (length(venn_up_membership_list) > 0) {
  venn_up_membership <- dplyr::bind_rows(
    venn_up_membership_list
  )
} else {
  venn_up_membership <- tibble::tibble(
    biotype_priority = character(),
    gene_id = character(),
    gene_name_priority = character(),
    HB = logical(),
    BT = logical(),
    HT = logical()
  )
}

write.table(
  venn_up_membership,
  file = file.path(
    venn_up_dir,
    "Venn_up_genes_membership_by_biotype.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Diagrammes de Venn terminés.")
message("Résultats dans : ", venn_up_dir)


## ----figures-19-venn-up-genes-by-biotype, echo=FALSE, message=FALSE, warning=FALSE, results="asis"-----------------------------------
fig_dir <- file.path(
  outdir,
  "Venn_gene_level_SPZ_up_by_biotype"
)

figs <- if (dir.exists(fig_dir)) {
  list.files(
    fig_dir,
    pattern = "^Venn_up_genes_.*\\.pdf$",
    full.names = TRUE
  )
} else {
  character(0)
}

include_pdf_as_large_png(figs)


## ----20-schematic-piecharts-by-biotype, eval=TRUE, message=FALSE, warning=FALSE------------------------------------------------------
############################################################
# 20. Schematic transition figure with pie charts by biotype
############################################################

message("---- Figures schématiques avec pie charts par biotype ----")

## The ggforce package is installed in block 1.
if (!requireNamespace("ggforce", quietly = TRUE)) {
  stop("Le package ggforce n'est pas disponible après le bloc d'installation.")
}

schematic_dir <- file.path(
  outdir,
  "Schematic_piecharts_by_biotype"
)

dir.create(
  schematic_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

## Remove old figures
old_schematic_files <- list.files(
  schematic_dir,
  pattern = "^Schematic_piecharts_.*\\.pdf$",
  full.names = TRUE
)

if (length(old_schematic_files) > 0) {
  invisible(file.remove(old_schematic_files))
}

preferred_biotype_order_schematic <- c(
  "miR",
  "tRF",
  "rRF",
  "snoRNA",
  "piR",
  "ncRNA",
  "cdna"
)

biotype_labels_pretty <- c(
  miR = "miRNA",
  tRF = "tRF",
  rRF = "rRF",
  snoRNA = "snoRNA",
  piR = "piRNA",
  ncRNA = "ncRNA",
  cdna = "cDNA"
)

contrast_labels_schematic <- c(
  Head_vs_Body = "H → B",
  Body_vs_Tail = "B → T",
  Head_vs_Tail = "H → T"
)

missing_contrasts_schematic <- setdiff(
  names(contrast_labels_schematic),
  names(gene_level_volcano_results)
)

if (length(missing_contrasts_schematic) > 0) {
  stop(
    "Contrastes manquants dans gene_level_volcano_results : ",
    paste(missing_contrasts_schematic, collapse = ", ")
  )
}

## Summarize DE genes and upregulated DE genes for each biotype
schematic_summary_raw <- lapply(
  names(contrast_labels_schematic),
  function(contrast_name) {

    df <- gene_level_volcano_results[[contrast_name]]

    if (is.null(df)) {
      stop("Résultats absents pour le contraste : ", contrast_name)
    }

    df %>%
      dplyr::filter(
        !is.na(biotype_priority),
        biotype_priority != "",
        !is.na(gene_id),
        gene_id != ""
      ) %>%
      dplyr::group_by(biotype_priority) %>%
      dplyr::summarise(
        n_DE_total = dplyr::n_distinct(
          gene_id[signif %in% TRUE]
        ),
        n_up = dplyr::n_distinct(
          gene_id[
            signif %in% TRUE &
              !is.na(log2FoldChange) &
              log2FoldChange < 0
          ]
        ),
        n_down = dplyr::n_distinct(
          gene_id[
            signif %in% TRUE &
              !is.na(log2FoldChange) &
              log2FoldChange > 0
          ]
        ),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        contrast = contrast_name,
        comparison = unname(
          contrast_labels_schematic[contrast_name]
        ),
        biotype_priority = as.character(biotype_priority)
      )
  }
) %>%
  dplyr::bind_rows()

observed_biotypes_schematic <- annot_gene %>%
  dplyr::filter(
    !is.na(biotype_priority),
    biotype_priority != ""
  ) %>%
  dplyr::distinct(biotype_priority) %>%
  dplyr::pull(biotype_priority) %>%
  as.character()

biotypes_schematic <- c(
  intersect(
    preferred_biotype_order_schematic,
    observed_biotypes_schematic
  ),
  setdiff(
    sort(observed_biotypes_schematic),
    preferred_biotype_order_schematic
  )
)

schematic_summary <- tidyr::expand_grid(
  biotype_priority = biotypes_schematic,
  contrast = names(contrast_labels_schematic)
) %>%
  dplyr::left_join(
    schematic_summary_raw,
    by = c("biotype_priority", "contrast")
  ) %>%
  dplyr::mutate(
    comparison = dplyr::coalesce(
      comparison,
      unname(contrast_labels_schematic[contrast])
    ),
    n_DE_total = tidyr::replace_na(n_DE_total, 0L),
    n_up = tidyr::replace_na(n_up, 0L),
    n_down = tidyr::replace_na(n_down, 0L)
  )

write.table(
  schematic_summary,
  file = file.path(
    schematic_dir,
    "Schematic_piecharts_counts_by_biotype.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## Positions of the three pie charts
pie_positions <- tibble::tibble(
  comparison = c("H → B", "B → T", "H → T"),
  x = c(2, 6, 4),
  y = c(1.15, 1.15, -1.55),
  r = c(0.55, 0.55, 0.62)
)

build_pie_df <- function(df_positioned) {

  pie_list <- lapply(
    seq_len(nrow(df_positioned)),
    function(i) {

      row_i <- df_positioned[i, , drop = FALSE]

      total_i <- as.numeric(row_i$n_DE_total)
      up_i <- as.numeric(row_i$n_up)
      other_i <- total_i - up_i

      ## No DE genes: gray circle with 0/0
      if (total_i == 0) {
        return(
          tibble::tibble(
            comparison = row_i$comparison,
            x = row_i$x,
            y = row_i$y,
            r = row_i$r,
            start = 0,
            end = 2 * pi,
            slice = "No DE genes"
          )
        )
      }

      slices_i <- list()

      if (up_i > 0) {
        slices_i[[length(slices_i) + 1]] <- tibble::tibble(
          comparison = row_i$comparison,
          x = row_i$x,
          y = row_i$y,
          r = row_i$r,
          start = 0,
          end = 2 * pi * up_i / total_i,
          slice = "Up genes"
        )
      }

      if (other_i > 0) {
        slices_i[[length(slices_i) + 1]] <- tibble::tibble(
          comparison = row_i$comparison,
          x = row_i$x,
          y = row_i$y,
          r = row_i$r,
          start = 2 * pi * up_i / total_i,
          end = 2 * pi,
          slice = "Other DE genes"
        )
      }

      dplyr::bind_rows(slices_i)
    }
  )

  dplyr::bind_rows(pie_list)
}

plot_one_biotype_schematic <- function(biotype_name) {

  df_bio <- schematic_summary %>%
    dplyr::filter(biotype_priority == biotype_name) %>%
    dplyr::left_join(
      pie_positions,
      by = "comparison"
    ) %>%
    dplyr::mutate(
      pie_title_y = dplyr::case_when(
        comparison == "H → T" ~ y + r + 0.28,
        TRUE ~ y + r + 0.22
      ),
      count_label_y = dplyr::case_when(
        comparison == "H → T" ~ y - r - 0.35,
        TRUE ~ y - r - 0.48
      )
    )

  pie_df <- build_pie_df(df_bio)

  if (biotype_name %in% names(biotype_labels_pretty)) {
    pretty_biotype <- unname(
      biotype_labels_pretty[biotype_name]
    )
  } else {
    pretty_biotype <- biotype_name
  }

  p <- ggplot() +

    ## Head to Body arrow
    annotate(
      "segment",
      x = 0.6,
      xend = 3.4,
      y = 0,
      yend = 0,
      arrow = grid::arrow(
        length = grid::unit(0.22, "cm"),
        type = "closed"
      ),
      linewidth = 0.8
    ) +

    ## Body to Tail arrow
    annotate(
      "segment",
      x = 4.6,
      xend = 7.4,
      y = 0,
      yend = 0,
      arrow = grid::arrow(
        length = grid::unit(0.22, "cm"),
        type = "closed"
      ),
      linewidth = 0.8
    ) +

    ## Curved Head to Tail arrow positioned at the bottom
    annotate(
      "curve",
      x = 0.7,
      y = -0.35,
      xend = 7.3,
      yend = -0.35,
      curvature = 0.45,
      arrow = grid::arrow(
        length = grid::unit(0.22, "cm"),
        type = "closed"
      ),
      linewidth = 0.8
    ) +

    ## Head, Body, and Tail nodes
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "H",
      size = 7,
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = 4,
      y = 0,
      label = "B",
      size = 7,
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = 8,
      y = 0,
      label = "T",
      size = 7,
      fontface = "bold"
    ) +

    ## Pie charts positioned near each arrow
    ggforce::geom_arc_bar(
      data = pie_df,
      aes(
        x0 = x,
        y0 = y,
        r0 = 0,
        r = r,
        start = start,
        end = end,
        fill = slice
      ),
      color = "black",
      linewidth = 0.35,
      inherit.aes = FALSE
    ) +

    ## Comparison name above the pie chart
    geom_text(
      data = df_bio,
      aes(
        x = x,
        y = pie_title_y,
        label = comparison
      ),
      fontface = "bold",
      size = 4.2
    ) +

    ## Up count / total DE count, separated from lines and arrows
    geom_label(
      data = df_bio,
      aes(
        x = x,
        y = count_label_y,
        label = paste0(n_up, "/", n_DE_total)
      ),
      size = 3.6,
      label.size = 0.2,
      fill = "white",
      alpha = 0.9
    ) +

    scale_fill_manual(
      values = c(
        "Up genes" = "black",
        "Other DE genes" = "grey82",
        "No DE genes" = "white"
      ),
      breaks = c(
        "Up genes",
        "Other DE genes",
        "No DE genes"
      ),
      labels = c(
        "Up genes",
        "Other DE genes",
        "No DE genes"
      ),
      drop = FALSE
    ) +

    coord_fixed(
      xlim = c(-0.8, 8.8),
      ylim = c(-2.8, 2.3),
      clip = "off"
    ) +

    labs(
      title = paste0(
        "Differential-expression transitions — ",
        pretty_biotype
      ),
      subtitle = paste0(
        "Black = genes higher in the second group (direction of the arrow); ",
        "labels = second-group-higher genes / all DE genes"
      ),
      fill = NULL
    ) +

    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      ),
      legend.position = "top"
    )

  biotype_file_name <- gsub(
    "[^A-Za-z0-9_-]+",
    "_",
    biotype_name
  )

  grDevices::pdf(
    file.path(
      schematic_dir,
      paste0(
        "Schematic_piecharts_",
        biotype_file_name,
        ".pdf"
      )
    ),
    width = 8,
    height = 5.6,
    useDingbats = FALSE
  )

  print(p)
  invisible(grDevices::dev.off())
}

for (bio in biotypes_schematic) {
  message("Figure schématique : ", bio)
  plot_one_biotype_schematic(bio)
}

message("Figures schématiques terminées.")
message("Résultats dans : ", schematic_dir)


## ----figures-20-schematic-piecharts-by-biotype, echo=FALSE, message=FALSE, warning=FALSE, results="asis"-----------------------------
fig_dir <- file.path(
  outdir,
  "Schematic_piecharts_by_biotype"
)

figs <- if (dir.exists(fig_dir)) {
  list.files(
    fig_dir,
    pattern = "^Schematic_piecharts_.*\\.pdf$",
    full.names = TRUE
  )
} else {
  character(0)
}

include_pdf_as_large_png(figs)


## ----21-relative-abundance-tail-enriched-de-small-rna, eval=TRUE, message=FALSE, warning=FALSE---------------------------------------
############################################################
# 21. Relative abundance of Tail-enriched DE small RNAs
############################################################

message(
  "---- Abondance relative des small RNA DE enrichis dans Tail ----"
)

relative_abundance_dir <- file.path(
  outdir,
  "Relative_abundance_Tail_enriched_DE"
)

dir.create(
  relative_abundance_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# Check required objects
############################################################

required_objects_abundance <- c(
  "counts_gene_DE",
  "annot_gene",
  "meta",
  "gene_level_volcano_results"
)

missing_objects_abundance <- required_objects_abundance[
  !vapply(
    required_objects_abundance,
    exists,
    logical(1)
  )
]

if (length(missing_objects_abundance) > 0) {
  stop(
    "Objets manquants pour le bloc d'abondance relative : ",
    paste(missing_objects_abundance, collapse = ", ")
  )
}

tail_contrasts <- c(
  "Head_vs_Tail",
  "Body_vs_Tail"
)

missing_tail_contrasts <- setdiff(
  tail_contrasts,
  names(gene_level_volcano_results)
)

if (length(missing_tail_contrasts) > 0) {
  stop(
    "Contrastes manquants dans gene_level_volcano_results : ",
    paste(missing_tail_contrasts, collapse = ", ")
  )
}

############################################################
# Identify Tail samples
############################################################

spz_tail_samples <- meta %>%
  dplyr::filter(
    `Tissue type` == "Spermatozoa",
    `Origin in Epididymis` == "Tail"
  ) %>%
  dplyr::pull(A_sample) %>%
  intersect(colnames(counts_gene_DE))

epididymis_tail_samples <- meta %>%
  dplyr::filter(
    `Tissue type` == "Epididymis",
    `Origin in Epididymis` == "Tail"
  ) %>%
  dplyr::pull(A_sample) %>%
  intersect(colnames(counts_gene_DE))

message(
  "Nombre d'échantillons spermatozoïdes Tail : ",
  length(spz_tail_samples)
)

message(
  "Nombre d'échantillons épididyme Tail : ",
  length(epididymis_tail_samples)
)

if (length(spz_tail_samples) == 0) {
  stop("Aucun échantillon spermatozoïde Tail trouvé.")
}

if (length(epididymis_tail_samples) == 0) {
  stop("Aucun échantillon épididyme Tail trouvé.")
}

############################################################
# Calculate gene-level relative abundance
############################################################

## Total annotated gene-level UMIs in each sample
gene_level_library_sizes <- colSums(
  counts_gene_DE,
  na.rm = TRUE
)

if (any(gene_level_library_sizes <= 0)) {
  stop(
    "Au moins un échantillon possède une taille de librairie ",
    "gene-level nulle."
  )
}

## Convert to UMI CPM
counts_gene_cpm <- sweep(
  counts_gene_DE,
  MARGIN = 2,
  STATS = gene_level_library_sizes,
  FUN = "/"
) * 1e6

############################################################
# Mean and median abundance in the two Tail sample types
############################################################

abundance_tail <- tibble::tibble(
  gene_id = rownames(counts_gene_DE),

  ## Tail spermatozoa: raw UMI counts
  mean_spz_tail_raw_umi = rowMeans(
    counts_gene_DE[
      ,
      spz_tail_samples,
      drop = FALSE
    ],
    na.rm = TRUE
  ),

  median_spz_tail_raw_umi = apply(
    counts_gene_DE[
      ,
      spz_tail_samples,
      drop = FALSE
    ],
    1,
    stats::median,
    na.rm = TRUE
  ),

  ## Tail spermatozoa: UMI CPM
  mean_spz_tail_cpm = rowMeans(
    counts_gene_cpm[
      ,
      spz_tail_samples,
      drop = FALSE
    ],
    na.rm = TRUE
  ),

  median_spz_tail_cpm = apply(
    counts_gene_cpm[
      ,
      spz_tail_samples,
      drop = FALSE
    ],
    1,
    stats::median,
    na.rm = TRUE
  ),

  ## Tail epididymis: raw UMI counts
  mean_epididymis_tail_raw_umi = rowMeans(
    counts_gene_DE[
      ,
      epididymis_tail_samples,
      drop = FALSE
    ],
    na.rm = TRUE
  ),

  median_epididymis_tail_raw_umi = apply(
    counts_gene_DE[
      ,
      epididymis_tail_samples,
      drop = FALSE
    ],
    1,
    stats::median,
    na.rm = TRUE
  ),

  ## Tail epididymis: UMI CPM
  mean_epididymis_tail_cpm = rowMeans(
    counts_gene_cpm[
      ,
      epididymis_tail_samples,
      drop = FALSE
    ],
    na.rm = TRUE
  ),

  median_epididymis_tail_cpm = apply(
    counts_gene_cpm[
      ,
      epididymis_tail_samples,
      drop = FALSE
    ],
    1,
    stats::median,
    na.rm = TRUE
  )
)

############################################################
# Abundance percentile within each biotype
############################################################
abundance_reference <- annot_gene %>%
  dplyr::filter(
    !is.na(gene_id),
    gene_id != "",
    !is.na(biotype_priority),
    biotype_priority != ""
  ) %>%
  dplyr::distinct(
    gene_id,
    .keep_all = TRUE
  ) %>%
  dplyr::left_join(
    abundance_tail,
    by = "gene_id"
  ) %>%
  dplyr::group_by(biotype_priority) %>%
  dplyr::mutate(

    ## A value close to 100 indicates a gene
    ## among the most abundant in its biotype.
    spz_tail_abundance_percentile = {
      if (dplyr::n() <= 1L) {
        rep(100, dplyr::n())
      } else {
        100 * dplyr::percent_rank(mean_spz_tail_cpm)
      }
    },

    epididymis_tail_abundance_percentile = {
      if (dplyr::n() <= 1L) {
        rep(100, dplyr::n())
      } else {
        100 * dplyr::percent_rank(
          mean_epididymis_tail_cpm
        )
      }
    }
  ) %>%
  dplyr::ungroup()

############################################################
# Orient DE results toward Tail
############################################################

comparison_labels_tail <- c(
  Head_vs_Tail = "Tail vs Head",
  Body_vs_Tail = "Tail vs Body"
)
tail_abundance_results <- lapply(
  tail_contrasts,
  function(contrast_name) {

    contrast_df <- gene_level_volcano_results[[contrast_name]]

    contrast_df %>%
      dplyr::filter(
        !is.na(gene_id),
        gene_id != "",
        !is.na(biotype_priority),
        biotype_priority != ""
      ) %>%
      dplyr::mutate(

        comparison = unname(
          comparison_labels_tail[contrast_name]
        ),

        ## The original contrasts are:
        ## Head vs Tail and Body vs Tail.
        ## Inverting the sign yields:
        ## Tail vs Head and Tail vs Body.
        tail_log2FoldChange = -log2FoldChange,

        tail_up_DE =
          signif %in% TRUE &
          !is.na(tail_log2FoldChange) &
          tail_log2FoldChange > 0,

        tail_down_DE =
          signif %in% TRUE &
          !is.na(tail_log2FoldChange) &
          tail_log2FoldChange < 0,

        DE_status = dplyr::case_when(
          tail_up_DE ~ "Tail-up DE",
          tail_down_DE ~ "Tail-down DE",
          TRUE ~ "Not DE"
        )
      ) %>%
      dplyr::left_join(
        abundance_reference %>%
          dplyr::select(
            gene_id,
            mean_spz_tail_raw_umi,
            median_spz_tail_raw_umi,
            mean_spz_tail_cpm,
            median_spz_tail_cpm,
            mean_epididymis_tail_raw_umi,
            median_epididymis_tail_raw_umi,
            mean_epididymis_tail_cpm,
            median_epididymis_tail_cpm,
            spz_tail_abundance_percentile,
            epididymis_tail_abundance_percentile
          ),
        by = "gene_id"
      )
  }
) %>%
  dplyr::bind_rows()

############################################################
# Order of comparisons and biotypes
############################################################

tail_abundance_results <- tail_abundance_results %>%
  dplyr::mutate(
    comparison = factor(
      comparison,
      levels = c(
        "Tail vs Head",
        "Tail vs Body"
      )
    )
  )

preferred_biotype_order_abundance <- c(
  "miR",
  "tRF",
  "rRF",
  "snoRNA",
  "piR",
  "ncRNA",
  "cdna"
)

observed_biotypes_abundance <- tail_abundance_results %>%
  dplyr::filter(tail_up_DE) %>%
  dplyr::distinct(biotype_priority) %>%
  dplyr::pull(biotype_priority) %>%
  as.character()

biotype_order_abundance <- c(
  intersect(
    preferred_biotype_order_abundance,
    observed_biotypes_abundance
  ),
  setdiff(
    sort(observed_biotypes_abundance),
    preferred_biotype_order_abundance
  )
)

if (length(biotype_order_abundance) == 0) {
  stop(
    "Aucun small RNA DE enrichi dans Tail n'a été trouvé ",
    "avec les seuils actuels."
  )
}

tail_abundance_results <- tail_abundance_results %>%
  dplyr::filter(
    biotype_priority %in% biotype_order_abundance
  ) %>%
  dplyr::mutate(
    biotype_priority = factor(
      biotype_priority,
      levels = biotype_order_abundance
    )
  )

############################################################
# Export the complete table
############################################################

write.table(
  tail_abundance_results,
  file = file.path(
    relative_abundance_dir,
    "Tail_enrichment_and_relative_abundance_all_genes.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Table of Tail-up DE genes
############################################################

tail_up_DE_abundance <- tail_abundance_results %>%
  dplyr::filter(tail_up_DE) %>%
  dplyr::arrange(
    biotype_priority,
    comparison,
    dplyr::desc(mean_spz_tail_cpm)
  )

write.table(
  tail_up_DE_abundance,
  file = file.path(
    relative_abundance_dir,
    "Tail_up_DE_genes_relative_abundance.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Abundance summary by biotype
############################################################

tail_up_abundance_summary <- tail_up_DE_abundance %>%
  dplyr::group_by(
    comparison,
    biotype_priority
  ) %>%
  dplyr::summarise(

    n_tail_up_DE = dplyr::n_distinct(gene_id),

    median_spz_tail_cpm = stats::median(
      mean_spz_tail_cpm,
      na.rm = TRUE
    ),

    mean_spz_tail_cpm = mean(
      mean_spz_tail_cpm,
      na.rm = TRUE
    ),

    max_spz_tail_cpm = max(
      mean_spz_tail_cpm,
      na.rm = TRUE
    ),

    n_in_top_10_percent_spz_tail = sum(
      spz_tail_abundance_percentile >= 90,
      na.rm = TRUE
    ),

    median_epididymis_tail_cpm = stats::median(
      mean_epididymis_tail_cpm,
      na.rm = TRUE
    ),

    mean_epididymis_tail_cpm = mean(
      mean_epididymis_tail_cpm,
      na.rm = TRUE
    ),

    max_epididymis_tail_cpm = max(
      mean_epididymis_tail_cpm,
      na.rm = TRUE
    ),

    n_in_top_10_percent_epididymis_tail = sum(
      epididymis_tail_abundance_percentile >= 90,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

print(tail_up_abundance_summary)

write.table(
  tail_up_abundance_summary,
  file = file.path(
    relative_abundance_dir,
    "Tail_up_DE_relative_abundance_summary_by_biotype.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Genes to label in the figures
############################################################

## Three most abundant Tail-up genes
## per biotype and comparison.
top_labels_spz_tail <- tail_up_DE_abundance %>%
  dplyr::filter(
    !is.na(gene_name_priority),
    gene_name_priority != ""
  ) %>%
  dplyr::group_by(
    comparison,
    biotype_priority
  ) %>%
  dplyr::slice_max(
    order_by = mean_spz_tail_cpm,
    n = 3,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

top_labels_epididymis_tail <- tail_up_DE_abundance %>%
  dplyr::filter(
    !is.na(gene_name_priority),
    gene_name_priority != ""
  ) %>%
  dplyr::group_by(
    comparison,
    biotype_priority
  ) %>%
  dplyr::slice_max(
    order_by = mean_epididymis_tail_cpm,
    n = 3,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

biotype_labels_abundance <- c(
  miR = "miRNA",
  tRF = "tRF",
  rRF = "rRF",
  snoRNA = "snoRNA",
  piR = "piRNA",
  ncRNA = "ncRNA",
  cdna = "cDNA"
)

############################################################
# Function for creating MA-like plots
############################################################

make_tail_abundance_MA_like_plot <- function(
  plot_df,
  label_df,
  abundance_column,
  x_axis_label,
  plot_title,
  output_filename
) {

  n_biotypes_plot <- dplyr::n_distinct(
    plot_df$biotype_priority
  )

  grDevices::pdf(
    file.path(
      relative_abundance_dir,
      output_filename
    ),
    width = 13,
    height = max(
      9,
      2.2 * n_biotypes_plot + 2
    ),
    useDingbats = FALSE
  )

  p <- ggplot(
    plot_df,
    aes(
      x = log10(.data[[abundance_column]] + 1),
      y = tail_log2FoldChange
    )
  ) +

    ## All genes in the biotype
    geom_point(
      color = "grey75",
      alpha = 0.45,
      size = 1
    ) +

    ## DE genes enriched in Tail
    geom_point(
      data = plot_df %>%
        dplyr::filter(tail_up_DE),
      color = "black",
      alpha = 0.90,
      size = 1.8
    ) +

    ## Line corresponding to no difference
    geom_hline(
      yintercept = 0,
      linewidth = 0.35
    ) +

    ## Log2 fold-change thresholds
    geom_hline(
      yintercept = c(
        -lfc_cut,
        lfc_cut
      ),
      linetype = "dashed",
      linewidth = 0.45
    ) +

    ## Label the most abundant Tail-up genes
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = gene_name_priority),
      size = 2.8,
      max.overlaps = 100,
      box.padding = 0.30,
      point.padding = 0.20,
      min.segment.length = 0
    ) +

    facet_grid(
      rows = ggplot2::vars(biotype_priority),
      cols = ggplot2::vars(comparison),
      scales = "free",
      labeller = ggplot2::labeller(
        biotype_priority = ggplot2::as_labeller(
          biotype_labels_abundance,
          default = ggplot2::label_value
        )
      )
    ) +

    labs(
      title = plot_title,
      subtitle = paste0(
        "Black = Tail-up DE genes | ",
        "positive log2FC = higher in Tail | ",
        "padj < ",
        padj_cut,
        " and |log2FC| ≥ ",
        lfc_cut
      ),
      x = x_axis_label,
      y = "log2 fold change: Tail vs comparison"
    ) +

    theme_bw() +

    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      ),
      strip.text = element_text(
        face = "bold"
      ),
      strip.text.y = element_text(
        angle = 0
      ),
      panel.grid.minor = element_blank()
    )

  print(p)
  invisible(grDevices::dev.off())
}

############################################################
# Figure 1: abundance in Tail spermatozoa
############################################################

make_tail_abundance_MA_like_plot(
  plot_df = tail_abundance_results,
  label_df = top_labels_spz_tail,
  abundance_column = "mean_spz_tail_cpm",
  x_axis_label = paste0(
    "log10(mean relative abundance in spermatozoa Tail ",
    "[UMI CPM] + 1)"
  ),
  plot_title = paste0(
    "Tail enrichment versus relative abundance ",
    "in spermatozoa Tail"
  ),
  output_filename = paste0(
    "MA_like_Tail_enrichment_",
    "SPZ_Tail_relative_abundance.pdf"
  )
)

############################################################
# Figure 2: abundance in Tail epididymis
############################################################

make_tail_abundance_MA_like_plot(
  plot_df = tail_abundance_results,
  label_df = top_labels_epididymis_tail,
  abundance_column = "mean_epididymis_tail_cpm",
  x_axis_label = paste0(
    "log10(mean relative abundance in epididymis Tail ",
    "[UMI CPM] + 1)"
  ),
  plot_title = paste0(
    "Tail enrichment in spermatozoa versus relative abundance ",
    "in epididymis Tail"
  ),
  output_filename = paste0(
    "MA_like_Tail_enrichment_",
    "Epididymis_Tail_relative_abundance.pdf"
  )
)

############################################################
# Figure 3: direct SPZ Tail / epididymis Tail comparison
############################################################

## The same gene can be Tail-up in one or both
## comparisons. Here, one row per gene is created.
tail_up_unique_genes <- tail_up_DE_abundance %>%
  dplyr::mutate(
    comparison_character = as.character(comparison)
  ) %>%
  dplyr::group_by(
    gene_id,
    gene_name_priority,
    biotype_priority
  ) %>%
  dplyr::summarise(
    comparisons_tail_up = paste(
      sort(unique(comparison_character)),
      collapse = "; "
    ),
    maximum_tail_log2FoldChange = max(
      tail_log2FoldChange,
      na.rm = TRUE
    ),
    mean_spz_tail_cpm = dplyr::first(
      mean_spz_tail_cpm
    ),
    mean_epididymis_tail_cpm = dplyr::first(
      mean_epididymis_tail_cpm
    ),
    .groups = "drop"
  )

top_labels_spz_vs_epididymis <- tail_up_unique_genes %>%
  dplyr::filter(
    !is.na(gene_name_priority),
    gene_name_priority != ""
  ) %>%
  dplyr::group_by(biotype_priority) %>%
  dplyr::slice_max(
    order_by = mean_spz_tail_cpm,
    n = 3,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

grDevices::pdf(
  file.path(
    relative_abundance_dir,
    paste0(
      "Scatter_SPZ_Tail_vs_Epididymis_Tail_",
      "abundance_of_Tail_up_DE_genes.pdf"
    )
  ),
  width = 12,
  height = 10,
  useDingbats = FALSE
)

p_spz_vs_epididymis_tail <- ggplot(
  tail_up_unique_genes,
  aes(
    x = log10(mean_spz_tail_cpm + 1),
    y = log10(mean_epididymis_tail_cpm + 1)
  )
) +

  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +

  geom_point(
    size = 2,
    alpha = 0.8
  ) +

  ggrepel::geom_text_repel(
    data = top_labels_spz_vs_epididymis,
    aes(label = gene_name_priority),
    size = 2.8,
    max.overlaps = 100,
    box.padding = 0.30,
    point.padding = 0.20,
    min.segment.length = 0
  ) +

  facet_wrap(
    ~ biotype_priority,
    scales = "free",
    labeller = ggplot2::as_labeller(
      biotype_labels_abundance,
      default = ggplot2::label_value
    )
  ) +

  labs(
    title = paste0(
      "Relative abundance of Tail-up DE small RNAs: ",
      "spermatozoa Tail versus epididymis Tail"
    ),
    subtitle = paste0(
      "Each point represents a unique gene enriched in ",
      "spermatozoa Tail"
    ),
    x = paste0(
      "log10(mean spermatozoa Tail abundance ",
      "[UMI CPM] + 1)"
    ),
    y = paste0(
      "log10(mean epididymis Tail abundance ",
      "[UMI CPM] + 1)"
    )
  ) +

  theme_bw() +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )

print(p_spz_vs_epididymis_tail)
invisible(grDevices::dev.off())

message("Analyse d'abondance relative terminée.")
message("Résultats dans : ", relative_abundance_dir)


## ----figures-21-relative-abundance-tail-enriched-de-small-rna, echo=FALSE, message=FALSE, warning=FALSE, results="asis"--------------
fig_dir <- file.path(
  outdir,
  "Relative_abundance_Tail_enriched_DE"
)

figs <- c(
  file.path(
    fig_dir,
    "MA_like_Tail_enrichment_SPZ_Tail_relative_abundance.pdf"
  ),
  file.path(
    fig_dir,
    "MA_like_Tail_enrichment_Epididymis_Tail_relative_abundance.pdf"
  ),
  file.path(
    fig_dir,
    paste0(
      "Scatter_SPZ_Tail_vs_Epididymis_Tail_",
      "abundance_of_Tail_up_DE_genes.pdf"
    )
  )
)

include_pdf_as_large_png(figs)


## ----22-read-level-gene-level-de-concordance-by-biotype, eval=TRUE, message=FALSE, warning=FALSE-------------------------------------
############################################################
# 22. Read-level and gene-level DE concordance by biotype
############################################################

message(
  "---- Concordance DE entre read-level et gene-level ----"
)

read_gene_status_dir <- file.path(
  outdir,
  "Read_gene_DE_status_by_biotype"
)

dir.create(
  read_gene_status_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# Check required objects and functions
############################################################

required_objects_read_gene_status <- c(
  "counts_DE",
  "counts_filtered_out",
  "meta",
  "gene_level_volcano_results",
  "padj_cut",
  "lfc_cut"
)

missing_objects_read_gene_status <-
  required_objects_read_gene_status[
    !vapply(
      required_objects_read_gene_status,
      exists,
      logical(1)
    )
  ]

if (length(missing_objects_read_gene_status) > 0) {
  stop(
    "Objets manquants pour l'analyse read/gene DE : ",
    paste(
      missing_objects_read_gene_status,
      collapse = ", "
    )
  )
}

if (!exists("prepare_spz_metadata", mode = "function")) {
  stop(
    "La fonction prepare_spz_metadata() est introuvable."
  )
}

if (!exists("clean_gene_id", mode = "function")) {
  stop(
    "La fonction clean_gene_id() est introuvable."
  )
}

############################################################
# Head / Body / Tail spermatozoa metadata
############################################################

meta_spz_read_status <- prepare_spz_metadata(
  regions = c(
    "Head",
    "Body",
    "Tail"
  )
)

meta_spz_read_status$Origin_in_Epididymis <- factor(
  meta_spz_read_status$Origin_in_Epididymis,
  levels = c(
    "Head",
    "Body",
    "Tail"
  )
)

counts_read_spz_status <- counts_DE[
  ,
  rownames(meta_spz_read_status),
  drop = FALSE
]

stopifnot(
  all(
    colnames(counts_read_spz_status) ==
      rownames(meta_spz_read_status)
  )
)

## Keep only read-level features present
## in at least one spermatozoa sample.
reads_present_in_spz <- rowSums(
  counts_read_spz_status,
  na.rm = TRUE
) > 0

counts_read_spz_status <- counts_read_spz_status[
  reads_present_in_spz,
  ,
  drop = FALSE
]

message(
  "Nombre de read-level features présentes dans les SPZ : ",
  nrow(counts_read_spz_status)
)

if (nrow(counts_read_spz_status) == 0) {
  stop(
    "Aucune read-level feature n'est présente dans les ",
    "échantillons spermatozoïdes."
  )
}

############################################################
# Shared read-level DESeq2 analysis across the three regions
############################################################

dds_read_spz_status <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(counts_read_spz_status),
  colData = meta_spz_read_status,
  design = ~ Origin_in_Epididymis
)

dds_read_spz_status <- DESeq2::DESeq(
  dds_read_spz_status
)

############################################################
# Read-level feature annotation
############################################################

read_annotation_status <- counts_filtered_out %>%
  dplyr::filter(
    feature_id %in% rownames(counts_read_spz_status),
    !is.na(gene_name_priority),
    gene_name_priority != "",
    gene_name_priority != "unannotated",
    !is.na(biotype_priority),
    biotype_priority != ""
  ) %>%
  dplyr::mutate(
    gene_id = clean_gene_id(
      gene_name_priority
    ),
    biotype_priority = as.character(
      biotype_priority
    )
  ) %>%
  dplyr::select(
    feature_id,
    sequence,
    gene_id,
    gene_name_priority,
    biotype_priority
  ) %>%
  dplyr::distinct(
    feature_id,
    .keep_all = TRUE
  )

message(
  "Nombre de read-level features annotées utilisées : ",
  nrow(read_annotation_status)
)

if (nrow(read_annotation_status) == 0) {
  stop(
    "Aucune read-level feature annotée n'est disponible."
  )
}

############################################################
# Define the three comparisons
############################################################

comparison_definitions_status <- list(
  Head_vs_Body = c(
    "Head",
    "Body"
  ),
  Body_vs_Tail = c(
    "Body",
    "Tail"
  ),
  Head_vs_Tail = c(
    "Head",
    "Tail"
  )
)

comparison_labels_status <- c(
  Head_vs_Body = "HB",
  Body_vs_Tail = "BT",
  Head_vs_Tail = "HT"
)

missing_gene_contrasts_status <- setdiff(
  names(comparison_definitions_status),
  names(gene_level_volcano_results)
)

if (length(missing_gene_contrasts_status) > 0) {
  stop(
    "Contrastes gene-level manquants : ",
    paste(
      missing_gene_contrasts_status,
      collapse = ", "
    )
  )
}

############################################################
# Compare read-level and gene-level statuses
############################################################

read_gene_status_results <- lapply(
  names(comparison_definitions_status),
  function(contrast_name) {

    groups_i <- comparison_definitions_status[[contrast_name]]

    group1 <- groups_i[1]
    group2 <- groups_i[2]

    message(
      "Analyse read/gene DE : ",
      contrast_name
    )

    ########################################################
    # Read-level results
    ########################################################

    read_res_i <- DESeq2::results(
      dds_read_spz_status,
      contrast = c(
        "Origin_in_Epididymis",
        group1,
        group2
      )
    )

    read_res_df_i <- as.data.frame(
      read_res_i
    ) %>%
      tibble::rownames_to_column(
        "feature_id"
      ) %>%
      dplyr::left_join(
        read_annotation_status,
        by = "feature_id"
      ) %>%
      dplyr::filter(
        !is.na(gene_id),
        gene_id != "",
        !is.na(biotype_priority),
        biotype_priority != ""
      ) %>%
      dplyr::mutate(
        read_DE =
          !is.na(padj) &
          padj < padj_cut &
          abs(log2FoldChange) >= lfc_cut
      ) %>%
      dplyr::rename(
        read_baseMean = baseMean,
        read_log2FoldChange = log2FoldChange,
        read_lfcSE = lfcSE,
        read_stat = stat,
        read_pvalue = pvalue,
        read_padj = padj
      )

    ########################################################
    # Gene-level results
    ########################################################

    gene_res_df_i <-
      gene_level_volcano_results[[contrast_name]] %>%
      dplyr::filter(
        !is.na(gene_id),
        gene_id != ""
      ) %>%
      dplyr::transmute(
        gene_id,
        gene_baseMean = baseMean,
        gene_log2FoldChange = log2FoldChange,
        gene_lfcSE = lfcSE,
        gene_stat = stat,
        gene_pvalue = pvalue,
        gene_padj = padj,
        gene_DE =
          !is.na(padj) &
          padj < padj_cut &
          abs(log2FoldChange) >= lfc_cut
      ) %>%
      dplyr::distinct(
        gene_id,
        .keep_all = TRUE
      )

    ########################################################
    # Associate each read with its gene
    ########################################################

    status_df_i <- read_res_df_i %>%
      dplyr::left_join(
        gene_res_df_i,
        by = "gene_id"
      )

    n_missing_gene_status <- sum(
      is.na(status_df_i$gene_DE)
    )

    if (n_missing_gene_status > 0) {
      stop(
        n_missing_gene_status,
        " read-level features n'ont pas trouvé leur ",
        "gène dans les résultats gene-level pour ",
        contrast_name,
        "."
      )
    }

    status_df_i %>%
      dplyr::mutate(
        contrast = contrast_name,

        comparison = unname(
          comparison_labels_status[contrast_name]
        ),

        status_category = dplyr::case_when(
          read_DE & gene_DE ~
            "Read DE / Gene DE",

          read_DE & !gene_DE ~
            "Read DE / Gene non-DE",

          !read_DE & gene_DE ~
            "Read non-DE / Gene DE",

          !read_DE & !gene_DE ~
            "Neither DE",

          TRUE ~ NA_character_
        )
      )
  }
) %>%
  dplyr::bind_rows()

############################################################
# Checks
############################################################

if (nrow(read_gene_status_results) == 0) {
  stop(
    "Le tableau de concordance read/gene est vide."
  )
}

if (any(is.na(read_gene_status_results$status_category))) {
  stop(
    "Au moins une read-level feature n'a pas reçu ",
    "de catégorie DE."
  )
}

message(
  "Nombre total de lignes read x contraste : ",
  nrow(read_gene_status_results)
)

############################################################
# Export the status of each read-level feature
############################################################

write.table(
  read_gene_status_results,
  file = file.path(
    read_gene_status_dir,
    paste0(
      "Read_level_and_gene_level_DE_status_",
      "all_features.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Number of DE reads out of the total number of reads per gene
############################################################

per_gene_read_DE_summary <- read_gene_status_results %>%
  dplyr::group_by(
    contrast,
    comparison,
    gene_id,
    gene_name_priority,
    biotype_priority,
    gene_DE
  ) %>%
  dplyr::summarise(

    n_total_reads = dplyr::n_distinct(
      feature_id
    ),

    n_DE_reads = dplyr::n_distinct(
      feature_id[
        read_DE %in% TRUE
      ]
    ),

    .groups = "drop"
  ) %>%
  dplyr::mutate(
    n_non_DE_reads =
      n_total_reads - n_DE_reads,

    percent_DE_reads = dplyr::if_else(
      n_total_reads > 0,
      100 * n_DE_reads / n_total_reads,
      0
    )
  ) %>%
  dplyr::arrange(
    contrast,
    biotype_priority,
    dplyr::desc(percent_DE_reads),
    dplyr::desc(n_total_reads)
  )

print(
  dplyr::slice_head(
    per_gene_read_DE_summary,
    n = 20
  )
)

write.table(
  per_gene_read_DE_summary,
  file = file.path(
    read_gene_status_dir,
    paste0(
      "Number_of_DE_reads_over_total_reads_",
      "within_each_gene.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Category summary by biotype
############################################################

status_levels_read_gene <- c(
  "Read DE / Gene DE",
  "Read DE / Gene non-DE",
  "Read non-DE / Gene DE",
  "Neither DE"
)

preferred_biotype_order_read_gene <- c(
  "miR",
  "tRF",
  "rRF",
  "snoRNA",
  "piR",
  "ncRNA",
  "cdna"
)

observed_biotypes_read_gene <-
  read_gene_status_results %>%
  dplyr::distinct(
    biotype_priority
  ) %>%
  dplyr::pull(
    biotype_priority
  ) %>%
  as.character()

biotype_order_read_gene <- c(
  intersect(
    preferred_biotype_order_read_gene,
    observed_biotypes_read_gene
  ),
  setdiff(
    sort(observed_biotypes_read_gene),
    preferred_biotype_order_read_gene
  )
)

read_gene_status_counts_raw <-
  read_gene_status_results %>%
  dplyr::count(
    contrast,
    comparison,
    biotype_priority,
    status_category,
    name = "n_reads"
  )

read_gene_status_barplot_summary <-
  tidyr::expand_grid(
    contrast = names(
      comparison_definitions_status
    ),
    biotype_priority =
      biotype_order_read_gene,
    status_category =
      status_levels_read_gene
  ) %>%
  dplyr::mutate(
    comparison = unname(
      comparison_labels_status[contrast]
    )
  ) %>%
  dplyr::left_join(
    read_gene_status_counts_raw,
    by = c(
      "contrast",
      "comparison",
      "biotype_priority",
      "status_category"
    )
  ) %>%
  dplyr::mutate(
    n_reads = tidyr::replace_na(
      n_reads,
      0L
    ),

    status_category = factor(
      status_category,
      levels = status_levels_read_gene
    ),

    biotype_priority = factor(
      biotype_priority,
      levels = biotype_order_read_gene
    ),

    comparison = factor(
      comparison,
      levels = c(
        "HB",
        "BT",
        "HT"
      )
    )
  ) %>%
  dplyr::group_by(
    contrast,
    comparison,
    biotype_priority
  ) %>%
  dplyr::mutate(
    n_total_reads_biotype = sum(
      n_reads
    ),

    percent_reads = dplyr::if_else(
      n_total_reads_biotype > 0,
      100 * n_reads / n_total_reads_biotype,
      0
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    comparison,
    biotype_priority,
    status_category
  )

print(read_gene_status_barplot_summary)

write.table(
  read_gene_status_barplot_summary,
  file = file.path(
    read_gene_status_dir,
    paste0(
      "Read_gene_DE_status_summary_",
      "by_biotype_and_comparison.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Names used in the figure
############################################################

biotype_labels_read_gene <- c(
  miR = "miRNA",
  tRF = "tRF",
  rRF = "rRF",
  snoRNA = "snoRNA",
  piR = "piRNA",
  ncRNA = "ncRNA",
  cdna = "cDNA"
)

status_labels_read_gene <- c(
  "Read DE / Gene DE" = "Read DE / Gene DE",
  "Read DE / Gene non-DE" = "Read DE / Gene non-DE",
  "Read non-DE / Gene DE" = "Read non-DE / Gene DE",
  "Neither DE" = "Neither DE"
)

############################################################
# Muted colors and stacking order
############################################################

status_stack_order <- c(
  "Neither DE",
  "Read non-DE / Gene DE",
  "Read DE / Gene non-DE",
  "Read DE / Gene DE"
)

status_colors_read_gene <- c(
  "Neither DE" = "#E6E8EA",
  "Read non-DE / Gene DE" = "#B8C4D0",
  "Read DE / Gene non-DE" = "#7E93A8",
  "Read DE / Gene DE" = "#43576B"
)

############################################################
# Prepare data for the figure
############################################################

plot_df_read_gene_stacked <- read_gene_status_barplot_summary %>%
  dplyr::mutate(
    biotype_label = dplyr::recode(
      as.character(biotype_priority),
      !!!biotype_labels_read_gene
    ),

    comparison = factor(
      as.character(comparison),
      levels = c("HB", "BT", "HT")
    ),

    status_category = factor(
      as.character(status_category),
      levels = status_stack_order
    )
  )

biotype_display_order <- c(
  "miRNA",
  "tRF",
  "rRF",
  "snoRNA",
  "piRNA",
  "ncRNA",
  "cDNA"
)

biotypes_present_in_plot <- biotype_display_order[
  biotype_display_order %in%
    unique(as.character(plot_df_read_gene_stacked$biotype_label))
]

plot_df_read_gene_stacked$biotype_label <- factor(
  plot_df_read_gene_stacked$biotype_label,
  levels = rev(biotypes_present_in_plot)
)

plot_df_read_gene_stacked <- plot_df_read_gene_stacked %>%
  dplyr::arrange(
    comparison,
    biotype_label,
    status_category
  ) %>%
  dplyr::group_by(
    comparison,
    biotype_label
  ) %>%
  dplyr::mutate(
    xmin = cumsum(dplyr::lag(percent_reads, default = 0)),
    xmax = xmin + percent_reads,
    xmid = xmin + percent_reads / 2,

    percent_label = dplyr::if_else(
      percent_reads >= 5,
      sprintf("%.1f%%", percent_reads),
      ""
    ),

    label_color = dplyr::if_else(
      status_category %in% c(
        "Read DE / Gene DE",
        "Read DE / Gene non-DE"
      ),
      "white",
      "black"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    y_center = as.numeric(biotype_label),
    ymin = y_center - 0.39,
    ymax = y_center + 0.39
  )

############################################################
# Single figure: one stacked bar per biotype
############################################################

output_file_stacked <- file.path(
  read_gene_status_dir,
  "Stacked_barplot_read_gene_DE_status_by_biotype.pdf"
)

grDevices::pdf(
  output_file_stacked,
  width = 15,
  height = 10,
  useDingbats = FALSE
)

p_read_gene_stacked <- ggplot(
  plot_df_read_gene_stacked
) +

  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = status_category
    ),
    color = "white",
    linewidth = 0.55
  ) +

  geom_text(
    data = plot_df_read_gene_stacked %>%
      dplyr::filter(percent_label != ""),
    aes(
      x = xmid,
      y = y_center,
      label = percent_label,
      color = label_color
    ),
    size = 3.2,
    fontface = "bold",
    show.legend = FALSE
  ) +

  facet_wrap(
    ~ comparison,
    ncol = 1
  ) +

  scale_fill_manual(
    values = status_colors_read_gene,
    breaks = status_stack_order,
    labels = status_labels_read_gene[status_stack_order],
    drop = FALSE
  ) +

  scale_color_identity() +

  scale_x_continuous(
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) {
      paste0(x, "%")
    },
    expand = c(0, 0)
  ) +
  
  coord_cartesian(
    xlim = c(0, 100)
  ) +

  scale_y_continuous(
    breaks = seq_along(levels(plot_df_read_gene_stacked$biotype_label)),
    labels = levels(plot_df_read_gene_stacked$biotype_label),
    expand = ggplot2::expansion(mult = c(0.06, 0.06))
  ) +

  labs(
    title = "Read-level and gene-level DE status by biotype",
    subtitle = paste0(
      "Each bar represents 100% of annotated read-level features ",
      "within one biotype | padj < ",
      padj_cut,
      " and |log2FC| ≥ ",
      lfc_cut
    ),
    x = "Percentage of read-level features",
    y = NULL,
    fill = NULL
  ) +

  theme_bw() +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    strip.background = element_rect(
      fill = "grey92",
      color = "grey50"
    ),
    legend.position = "top",
    legend.text = element_text(
      size = 9
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(
      face = "bold",
      size = 10
    ),
    axis.text.x = element_text(
      size = 9
    )
  )

print(p_read_gene_stacked)

invisible(grDevices::dev.off())

message(
  "Analyse read-level / gene-level terminée."
)

message(
  "Figure générée : ",
  output_file_stacked
)

message(
  "Résultats dans : ",
  read_gene_status_dir
)


## ----figures-22-read-level-gene-level-de-concordance-by-biotype, echo=FALSE, message=FALSE, warning=FALSE, results="asis"------------
figs <- c(
  file.path(
    outdir,
    "Read_gene_DE_status_by_biotype",
    "Stacked_barplot_read_gene_DE_status_by_biotype.pdf"
  )
)

include_pdf_as_large_png(figs)


## ----23-correlation-heatmaps-by-biotype, eval=TRUE, message=FALSE, warning=FALSE-----------------------------------------------------
############################################################
# 23. Sample correlation heatmaps by biotype
############################################################

message(
  "---- Heatmaps de corrélation par biotype ----"
)

############################################################
# Check the pheatmap package
############################################################

if (!requireNamespace("pheatmap", quietly = TRUE)) {

  cran_repository <- getOption("repos")[["CRAN"]]

  if (
    is.null(cran_repository) ||
    is.na(cran_repository) ||
    cran_repository == "@CRAN@"
  ) {
    cran_repository <- "https://cloud.r-project.org"
  }

  install.packages(
    "pheatmap",
    repos = cran_repository
  )
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop(
    "Le package pheatmap n'a pas pu être chargé."
  )
}

############################################################
# Output directory
############################################################

correlation_root_dir <- file.path(
  outdir,
  "Correlation_heatmaps_by_biotype"
)

dir.create(
  correlation_root_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# Check required objects
############################################################

required_objects_correlation <- c(
  "counts_DE",
  "counts_filtered_out",
  "meta"
)

missing_objects_correlation <- required_objects_correlation[
  !vapply(
    required_objects_correlation,
    exists,
    logical(1)
  )
]

if (length(missing_objects_correlation) > 0) {
  stop(
    "Objets manquants pour les heatmaps de corrélation : ",
    paste(
      missing_objects_correlation,
      collapse = ", "
    )
  )
}

if (!exists("safe_vst", mode = "function")) {
  stop(
    "La fonction safe_vst() est introuvable."
  )
}

############################################################
# Remove old figures
############################################################

old_correlation_files <- list.files(
  correlation_root_dir,
  pattern = "\\.(pdf|png)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(old_correlation_files) > 0) {
  invisible(
    file.remove(old_correlation_files)
  )
}

############################################################
# Biotype order and names
############################################################

preferred_biotype_order_correlation <- c(
  "miR",
  "tRF",
  "rRF",
  "snoRNA",
  "piR",
  "ncRNA",
  "cdna"
)

observed_biotypes_correlation <- counts_filtered_out %>%
  dplyr::filter(
    !is.na(biotype_priority),
    biotype_priority != "",
    !is.na(feature_id),
    feature_id != ""
  ) %>%
  dplyr::distinct(
    biotype_priority
  ) %>%
  dplyr::pull(
    biotype_priority
  ) %>%
  as.character()

biotypes_for_correlation <- c(
  intersect(
    preferred_biotype_order_correlation,
    observed_biotypes_correlation
  ),
  setdiff(
    sort(observed_biotypes_correlation),
    preferred_biotype_order_correlation
  )
)

if (length(biotypes_for_correlation) == 0) {
  stop(
    "Aucun biotype annoté n'est disponible pour les heatmaps."
  )
}

biotype_labels_correlation <- c(
  miR = "miRNA",
  tRF = "tRF",
  rRF = "rRF",
  snoRNA = "snoRNA",
  piR = "piRNA",
  ncRNA = "ncRNA",
  cdna = "cDNA"
)

############################################################
# Contexts analyzed separately
############################################################

correlation_contexts <- list(

  Spermatozoa = list(
    tissue = "Spermatozoa",
    folder = "Spermatozoa",
    file_prefix = "SPZ",
    title = "Spermatozoa"
  ),

  Epididymis = list(
    tissue = "Epididymis",
    folder = "Epididymis",
    file_prefix = "Epididymis",
    title = "Epididymis"
  )
)

############################################################
# Fixed palette from 0 to 1
############################################################

correlation_colors <- grDevices::colorRampPalette(
  c(
    "#2166AC",
    "#67A9CF",
    "white",
    "#EF8A62",
    "#B2182B"
  )
)(101)

## The number of breaks must equal
## the number of colors + 1.
correlation_breaks <- seq(
  from = 0,
  to = 1,
  length.out = length(correlation_colors) + 1
)

############################################################
# Function for creating a heatmap
############################################################

make_correlation_heatmap_one_biotype <- function(
  biotype_name,
  meta_context,
  samples_context,
  context_dir,
  context_file_prefix,
  context_title
) {

  message(
    "Heatmap de corrélation — ",
    context_title,
    " — biotype ",
    biotype_name
  )

  ##########################################################
  # Select features belonging to the biotype
  ##########################################################

  features_biotype <- counts_filtered_out %>%
    dplyr::filter(
      as.character(biotype_priority) == biotype_name,
      !is.na(feature_id),
      feature_id != ""
    ) %>%
    dplyr::pull(
      feature_id
    ) %>%
    unique() %>%
    intersect(
      rownames(counts_DE)
    )

  message(
    "Nombre initial de features : ",
    length(features_biotype)
  )

  if (length(features_biotype) < 2) {
    warning(
      "Biotype ",
      biotype_name,
      " ignoré dans ",
      context_title,
      " : moins de deux features."
    )

    return(
      invisible(NULL)
    )
  }

  ##########################################################
  # Matrix for the biotype and tissue of interest
  ##########################################################

  counts_biotype_context <- counts_DE[
    features_biotype,
    samples_context,
    drop = FALSE
  ]

  ## Remove features absent from all samples
  counts_biotype_context <- counts_biotype_context[
    rowSums(
      counts_biotype_context,
      na.rm = TRUE
    ) > 0,
    ,
    drop = FALSE
  ]

  message(
    "Nombre de features après filtre > 0 : ",
    nrow(counts_biotype_context)
  )

  if (nrow(counts_biotype_context) < 2) {
    warning(
      "Biotype ",
      biotype_name,
      " ignoré dans ",
      context_title,
      " : moins de deux features exprimées."
    )

    return(
      invisible(NULL)
    )
  }

  ## Check that no sample has an entirely zero column
  samples_all_zero <- colSums(
    counts_biotype_context,
    na.rm = TRUE
  ) == 0

  if (any(samples_all_zero)) {
    warning(
      "Biotype ",
      biotype_name,
      " ignoré dans ",
      context_title,
      " : aucun count dans les échantillons ",
      paste(
        colnames(counts_biotype_context)[samples_all_zero],
        collapse = ", "
      ),
      "."
    )

    return(
      invisible(NULL)
    )
  }

  ##########################################################
  # Strict metadata alignment
  ##########################################################

  meta_context_aligned <- meta_context[
    colnames(counts_biotype_context),
    ,
    drop = FALSE
  ]

  stopifnot(
    all(
      rownames(meta_context_aligned) ==
        colnames(counts_biotype_context)
    )
  )

  ##########################################################
  # VST transformation
  ##########################################################

  dds_correlation <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(
      counts_biotype_context
    ),
    colData = meta_context_aligned,
    design = ~ 1
  )

  ## More robust for matrices containing many zeros
  dds_correlation <- DESeq2::estimateSizeFactors(
    dds_correlation,
    type = "poscounts"
  )

  vsd_correlation <- safe_vst(
    dds_correlation,
    label = paste0(
      "correlation ",
      context_title,
      " — ",
      biotype_name
    )
  )

  vst_matrix_correlation <- SummarizedExperiment::assay(
    vsd_correlation
  )

  ##########################################################
  # Spearman correlation between samples
  ##########################################################

  correlation_matrix <- stats::cor(
    vst_matrix_correlation,
    method = "spearman",
    use = "pairwise.complete.obs"
  )

  if (any(!is.finite(correlation_matrix))) {
    warning(
      "Heatmap ignorée pour ",
      biotype_name,
      " dans ",
      context_title,
      " : la matrice contient des corrélations non finies."
    )

    return(
      invisible(NULL)
    )
  }

  ## Numerical safeguard
  correlation_matrix[
    correlation_matrix > 1
  ] <- 1

  correlation_matrix[
    correlation_matrix < -1
  ] <- -1

  ##########################################################
  # Sample labels
  ##########################################################

  meta_context_labels <- meta_context_aligned %>%
    as.data.frame() %>%
    dplyr::mutate(

      sample_base_label = dplyr::case_when(

        !is.na(`Unique Sample name`) &
          trimws(
            as.character(`Unique Sample name`)
          ) != "" ~

          as.character(`Unique Sample name`),

        TRUE ~
          as.character(A_sample)
      ),

      sample_display_label = paste0(
        sample_base_label,
        " | ",
        as.character(`Origin in Epididymis`)
      )
    )

  sample_label_vector <- stats::setNames(
    meta_context_labels$sample_display_label,
    meta_context_labels$A_sample
  )

  sample_ids <- colnames(
    correlation_matrix
  )

  display_labels <- unname(
    sample_label_vector[
      sample_ids
    ]
  )

  missing_labels <- is.na(
    display_labels
  )

  display_labels[
    missing_labels
  ] <- sample_ids[
    missing_labels
  ]

  ## Avoid duplicate names in the heatmap
  display_labels <- make.unique(
    display_labels
  )

  rownames(
    correlation_matrix
  ) <- display_labels

  colnames(
    correlation_matrix
  ) <- display_labels

  ##########################################################
  # Export the correlation matrix
  ##########################################################

  biotype_file_name <- gsub(
    "[^A-Za-z0-9_-]+",
    "_",
    biotype_name
  )

  correlation_table_out <- as.data.frame(
    correlation_matrix
  ) %>%
    tibble::rownames_to_column(
      "sample"
    )

  write.table(
    correlation_table_out,
    file = file.path(
      context_dir,
      paste0(
        "Correlation_matrix_",
        context_file_prefix,
        "_",
        biotype_file_name,
        ".tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  ##########################################################
  # Readable biotype name
  ##########################################################

  if (biotype_name %in% names(biotype_labels_correlation)) {

    pretty_biotype <- unname(
      biotype_labels_correlation[
        biotype_name
      ]
    )

  } else {

    pretty_biotype <- biotype_name
  }

  ##########################################################
  # Dimension parameters
  ##########################################################

  n_samples_correlation <- ncol(
    correlation_matrix
  )

  pdf_size_correlation <- max(
    9,
    0.50 * n_samples_correlation + 5
  )

  fontsize_samples <- dplyr::case_when(
    n_samples_correlation <= 12 ~ 8,
    n_samples_correlation <= 20 ~ 7,
    TRUE ~ 6
  )

  ##########################################################
  # Create the PDF directly with pheatmap
  ##########################################################

  output_pdf_correlation <- file.path(
    context_dir,
    paste0(
      "Correlation_heatmap_",
      context_file_prefix,
      "_",
      biotype_file_name,
      ".pdf"
    )
  )

  ## Remove any incomplete file
  if (file.exists(output_pdf_correlation)) {
    invisible(
      file.remove(output_pdf_correlation)
    )
  }

  heatmap_success <- tryCatch(
    {

      pheatmap::pheatmap(
        mat = correlation_matrix,

        color = correlation_colors,
        breaks = correlation_breaks,

        cluster_rows = TRUE,
        cluster_cols = TRUE,

        clustering_distance_rows = stats::as.dist(
          1 - correlation_matrix
        ),

        clustering_distance_cols = stats::as.dist(
          1 - correlation_matrix
        ),

        clustering_method = "average",

        ## No numbers in the cells
        display_numbers = FALSE,

        ## No visible borders between cells
        border_color = NA,

        ## Show dendrograms
        treeheight_row = 50,
        treeheight_col = 50,

        fontsize = 8,
        fontsize_row = fontsize_samples,
        fontsize_col = fontsize_samples,

        angle_col = 45,

        legend = TRUE,

        main = paste0(
          "Sample correlation — ",
          context_title,
          " — ",
          pretty_biotype
        ),

        filename = output_pdf_correlation,
        width = pdf_size_correlation,
        height = pdf_size_correlation,

        silent = TRUE
      )

      TRUE
    },

    error = function(e) {

      warning(
        "Échec de la heatmap ",
        context_title,
        " — ",
        biotype_name,
        " : ",
        conditionMessage(e)
      )

      FALSE
    }
  )

  ##########################################################
  # Check the generated PDF
  ##########################################################

  pdf_is_created <-
    heatmap_success &&
    file.exists(output_pdf_correlation) &&
    !is.na(file.info(output_pdf_correlation)$size) &&
    file.info(output_pdf_correlation)$size > 0

  if (!pdf_is_created) {

    if (file.exists(output_pdf_correlation)) {
      invisible(
        file.remove(output_pdf_correlation)
      )
    }

    warning(
      "Le PDF n'a pas été correctement généré : ",
      output_pdf_correlation
    )

    return(
      invisible(NULL)
    )
  }

  message(
    "Heatmap créée : ",
    output_pdf_correlation
  )

  invisible(
    output_pdf_correlation
  )
}

############################################################
# Generate heatmaps
############################################################

generated_correlation_files <- list()

for (context_name in names(correlation_contexts)) {

  context_info <- correlation_contexts[[context_name]]

  context_dir <- file.path(
    correlation_root_dir,
    context_info$folder
  )

  dir.create(
    context_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  ##########################################################
  # Metadata for the tissue of interest
  ##########################################################

  meta_context <- meta %>%
    dplyr::filter(
      `Tissue type` == context_info$tissue,
      `Origin in Epididymis` %in% c(
        "Head",
        "Body",
        "Tail"
      )
    )

  meta_context <- as.data.frame(
    meta_context
  )

  rownames(
    meta_context
  ) <- meta_context$A_sample

  samples_context <- intersect(
    rownames(meta_context),
    colnames(counts_DE)
  )

  meta_context <- meta_context[
    samples_context,
    ,
    drop = FALSE
  ]

  meta_context$`Origin in Epididymis` <- factor(
    meta_context$`Origin in Epididymis`,
    levels = c(
      "Head",
      "Body",
      "Tail"
    )
  )

  meta_context <- droplevels(
    meta_context
  )

  message(
    "Contexte ",
    context_name,
    " : ",
    length(samples_context),
    " échantillons."
  )

  print(
    table(
      meta_context$`Origin in Epididymis`,
      useNA = "ifany"
    )
  )

  if (length(samples_context) < 3) {
    warning(
      "Contexte ",
      context_name,
      " ignoré : moins de trois échantillons."
    )

    next
  }

  ##########################################################
  # One heatmap per biotype
  ##########################################################

  generated_correlation_files[[context_name]] <- lapply(
    biotypes_for_correlation,
    function(biotype_name) {

      make_correlation_heatmap_one_biotype(
        biotype_name = biotype_name,
        meta_context = meta_context,
        samples_context = samples_context,
        context_dir = context_dir,
        context_file_prefix = context_info$file_prefix,
        context_title = context_info$title
      )
    }
  )

  names(
    generated_correlation_files[[context_name]]
  ) <- biotypes_for_correlation
}

message(
  "Heatmaps de corrélation par biotype terminées."
)

message(
  "Résultats dans : ",
  correlation_root_dir
)


## ----figures-23-correlation-heatmaps-by-biotype, echo=FALSE, message=FALSE, warning=FALSE, results="asis"----------------------------
correlation_root_dir <- file.path(
  outdir,
  "Correlation_heatmaps_by_biotype"
)

figs_spz <- if (
  dir.exists(
    file.path(
      correlation_root_dir,
      "Spermatozoa"
    )
  )
) {

  list.files(
    file.path(
      correlation_root_dir,
      "Spermatozoa"
    ),
    pattern = "^Correlation_heatmap_SPZ_.*\\.pdf$",
    full.names = TRUE
  )

} else {

  character(0)
}

figs_epididymis <- if (
  dir.exists(
    file.path(
      correlation_root_dir,
      "Epididymis"
    )
  )
) {

  list.files(
    file.path(
      correlation_root_dir,
      "Epididymis"
    ),
    pattern = "^Correlation_heatmap_Epididymis_.*\\.pdf$",
    full.names = TRUE
  )

} else {

  character(0)
}

figs <- c(
  sort(figs_spz),
  sort(figs_epididymis)
)

include_pdf_as_large_png(
  figs,
  dpi = 180
)


## ----24-distribution-filtered-dataset-and-biotypes, eval=TRUE, message=FALSE, warning=FALSE------------------------------------------
############################################################
# 24. Distribution of the filtered dataset and biotypes
############################################################

message(
  "---- Distribution du dataset filtré et des biotypes ----"
)

distribution_dir <- file.path(
  outdir,
  "Distribution_filtered_dataset_and_biotypes"
)

dir.create(
  distribution_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# Check required objects
############################################################

required_objects_distribution <- c(
  "counts_DE",
  "counts_filtered_out",
  "meta"
)

missing_objects_distribution <- required_objects_distribution[
  !vapply(
    required_objects_distribution,
    exists,
    logical(1)
  )
]

if (length(missing_objects_distribution) > 0) {
  stop(
    "Objets manquants pour le bloc distribution : ",
    paste(
      missing_objects_distribution,
      collapse = ", "
    )
  )
}

if (!exists("safe_vst", mode = "function")) {
  stop(
    "La fonction safe_vst() est introuvable."
  )
}

############################################################
# Annotation of retained features
############################################################

distribution_annotation <- counts_filtered_out %>%
  dplyr::filter(
    feature_id %in% rownames(counts_DE),
    !is.na(feature_id),
    feature_id != "",
    !is.na(biotype_priority),
    biotype_priority != ""
  ) %>%
  dplyr::select(
    feature_id,
    biotype_priority
  ) %>%
  dplyr::distinct(
    feature_id,
    .keep_all = TRUE
  )

if (nrow(distribution_annotation) == 0) {
  stop(
    "Aucune feature annotée n'a été retrouvée dans counts_DE."
  )
}

############################################################
# VST on the entire filtered dataset
############################################################

meta_distribution <- meta %>%
  dplyr::filter(
    A_sample %in% colnames(counts_DE)
  )

meta_distribution <- as.data.frame(meta_distribution)
rownames(meta_distribution) <- meta_distribution$A_sample
meta_distribution <- meta_distribution[
  colnames(counts_DE),
  ,
  drop = FALSE
]

dds_distribution <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(counts_DE),
  colData = meta_distribution,
  design = ~ 1
)

dds_distribution <- DESeq2::estimateSizeFactors(
  dds_distribution,
  type = "poscounts"
)

vsd_distribution <- safe_vst(
  dds_distribution,
  label = "distribution filtered dataset"
)

vst_distribution_matrix <- SummarizedExperiment::assay(
  vsd_distribution
)

############################################################
# Summary table by feature
############################################################

feature_distribution_df <- tibble::tibble(
  feature_id = rownames(counts_DE),
  mean_raw_count = rowMeans(
    counts_DE,
    na.rm = TRUE
  ),
  log10_mean_raw_count = log10(
    rowMeans(
      counts_DE,
      na.rm = TRUE
    ) + 1
  ),
  mean_vst = rowMeans(
    vst_distribution_matrix,
    na.rm = TRUE
  )
) %>%
  dplyr::left_join(
    distribution_annotation,
    by = "feature_id"
  ) %>%
  dplyr::filter(
    !is.na(biotype_priority),
    biotype_priority != ""
  )

write.table(
  feature_distribution_df,
  file = file.path(
    distribution_dir,
    "Feature_distribution_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Biotype order
############################################################

preferred_biotype_order_distribution <- c(
  "miR",
  "tRF",
  "rRF",
  "snoRNA",
  "piR",
  "ncRNA",
  "cdna"
)

observed_biotypes_distribution <- feature_distribution_df %>%
  dplyr::distinct(biotype_priority) %>%
  dplyr::pull(biotype_priority) %>%
  as.character()

biotype_order_distribution <- c(
  intersect(
    preferred_biotype_order_distribution,
    observed_biotypes_distribution
  ),
  setdiff(
    sort(observed_biotypes_distribution),
    preferred_biotype_order_distribution
  )
)

biotype_labels_distribution <- c(
  miR = "miRNA",
  tRF = "tRF",
  rRF = "rRF",
  snoRNA = "snoRNA",
  piR = "piRNA",
  ncRNA = "ncRNA",
  cdna = "cDNA"
)

feature_distribution_df <- feature_distribution_df %>%
  dplyr::mutate(
    biotype_priority = factor(
      biotype_priority,
      levels = biotype_order_distribution
    )
  )

############################################################
# Figure 1: global histogram of log10(mean count + 1)
############################################################

grDevices::pdf(
  file.path(
    distribution_dir,
    "Histogram_global_log10_mean_raw_count.pdf"
  ),
  width = 10,
  height = 7,
  useDingbats = FALSE
)

p_hist_global <- ggplot(
  feature_distribution_df,
  aes(x = log10_mean_raw_count)
) +
  geom_histogram(
    bins = 60,
    fill = "grey70",
    color = "white"
  ) +
  geom_density(
    linewidth = 0.9
  ) +
  theme_bw() +
  labs(
    title = "Global distribution of filtered features",
    subtitle = "log10(mean raw count + 1)",
    x = "log10(mean raw count + 1)",
    y = "Number of features"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    panel.grid.minor = element_blank()
  )

print(p_hist_global)
invisible(grDevices::dev.off())

############################################################
# Figure 2: histograms of mean counts by biotype
############################################################

grDevices::pdf(
  file.path(
    distribution_dir,
    "Histogram_by_biotype_log10_mean_raw_count.pdf"
  ),
  width = 14,
  height = 10,
  useDingbats = FALSE
)

p_hist_biotype_raw <- ggplot(
  feature_distribution_df,
  aes(x = log10_mean_raw_count)
) +
  geom_histogram(
    bins = 40,
    fill = "grey70",
    color = "white"
  ) +
  facet_wrap(
    ~ biotype_priority,
    scales = "free_y",
    ncol = 3,
    labeller = ggplot2::as_labeller(
      biotype_labels_distribution,
      default = ggplot2::label_value
    )
  ) +
  theme_bw() +
  labs(
    title = "Distribution by biotype",
    subtitle = "log10(mean raw count + 1)",
    x = "log10(mean raw count + 1)",
    y = "Number of features"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )

print(p_hist_biotype_raw)
invisible(grDevices::dev.off())

############################################################
# Figure 3: histograms of mean VST by biotype
############################################################

grDevices::pdf(
  file.path(
    distribution_dir,
    "Histogram_by_biotype_mean_VST.pdf"
  ),
  width = 14,
  height = 10,
  useDingbats = FALSE
)

p_hist_biotype_vst <- ggplot(
  feature_distribution_df,
  aes(x = mean_vst)
) +
  geom_histogram(
    bins = 40,
    fill = "grey70",
    color = "white"
  ) +
  facet_wrap(
    ~ biotype_priority,
    scales = "free_y",
    ncol = 3,
    labeller = ggplot2::as_labeller(
      biotype_labels_distribution,
      default = ggplot2::label_value
    )
  ) +
  theme_bw() +
  labs(
    title = "Distribution by biotype",
    subtitle = "Mean VST per feature",
    x = "Mean VST",
    y = "Number of features"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )

print(p_hist_biotype_vst)
invisible(grDevices::dev.off())

############################################################
# Figure 4: QQ plots of mean VST by biotype
############################################################

grDevices::pdf(
  file.path(
    distribution_dir,
    "QQplot_by_biotype_mean_VST.pdf"
  ),
  width = 14,
  height = 10,
  useDingbats = FALSE
)

p_qq_biotype_vst <- ggplot(
  feature_distribution_df,
  aes(sample = mean_vst)
) +
  stat_qq(
    size = 0.6,
    alpha = 0.7
  ) +
  stat_qq_line(
    linewidth = 0.6
  ) +
  facet_wrap(
    ~ biotype_priority,
    scales = "free",
    ncol = 3,
    labeller = ggplot2::as_labeller(
      biotype_labels_distribution,
      default = ggplot2::label_value
    )
  ) +
  theme_bw() +
  labs(
    title = "QQ-plots by biotype",
    subtitle = "Mean VST per feature",
    x = "Theoretical quantiles",
    y = "Observed quantiles"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )

print(p_qq_biotype_vst)
invisible(grDevices::dev.off())

message(
  "Bloc distribution terminé."
)

message(
  "Résultats dans : ",
  distribution_dir
)


## ----figures-24-distribution-filtered-dataset-and-biotypes, echo=FALSE, message=FALSE, warning=FALSE, results="asis"-----------------
distribution_dir <- file.path(
  outdir,
  "Distribution_filtered_dataset_and_biotypes"
)

figs <- c(
  file.path(
    distribution_dir,
    "Histogram_global_log10_mean_raw_count.pdf"
  ),
  file.path(
    distribution_dir,
    "Histogram_by_biotype_log10_mean_raw_count.pdf"
  ),
  file.path(
    distribution_dir,
    "Histogram_by_biotype_mean_VST.pdf"
  ),
  file.path(
    distribution_dir,
    "QQplot_by_biotype_mean_VST.pdf"
  )
)

include_pdf_as_large_png(figs)


## ----22b-barplot-total-read-features-and-up-down-de-by-biotype, eval=TRUE, message=FALSE, warning=FALSE------------------------------
############################################################
# 22b. Total read-level features and up/down DE features
#         by biotype
############################################################

message(
  paste0(
    "---- Barplot du nombre total de read-level features ",
    "et des features DE up/down par biotype ----"
  )
)

read_barplot_de_dir <- file.path(
  outdir,
  "Barplot_DE_read_features_by_biotype"
)

dir.create(
  read_barplot_de_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# Check required objects
############################################################

required_objects_read_barplot <- c(
  "read_gene_status_results",
  "padj_cut",
  "lfc_cut"
)

missing_objects_read_barplot <- required_objects_read_barplot[
  !vapply(
    required_objects_read_barplot,
    exists,
    logical(1)
  )
]

if (length(missing_objects_read_barplot) > 0) {
  stop(
    "Objets manquants pour le barplot read-level : ",
    paste(
      missing_objects_read_barplot,
      collapse = ", "
    ),
    ". Le bloc 22 doit être exécuté avant ce bloc."
  )
}

required_columns_read_barplot <- c(
  "feature_id",
  "biotype_priority",
  "contrast",
  "comparison",
  "read_DE",
  "read_log2FoldChange"
)

missing_columns_read_barplot <- setdiff(
  required_columns_read_barplot,
  colnames(read_gene_status_results)
)

if (length(missing_columns_read_barplot) > 0) {
  stop(
    "Colonnes manquantes dans read_gene_status_results : ",
    paste(
      missing_columns_read_barplot,
      collapse = ", "
    )
  )
}

############################################################
# Comparison order
############################################################

contrast_order_read_barplot <- c(
  "Head_vs_Body",
  "Body_vs_Tail",
  "Head_vs_Tail"
)

comparison_labels_read_barplot <- c(
  Head_vs_Body = "HB",
  Body_vs_Tail = "BT",
  Head_vs_Tail = "HT"
)

missing_contrasts_read_barplot <- setdiff(
  contrast_order_read_barplot,
  unique(
    as.character(
      read_gene_status_results$contrast
    )
  )
)

if (length(missing_contrasts_read_barplot) > 0) {
  stop(
    "Contrastes read-level manquants : ",
    paste(
      missing_contrasts_read_barplot,
      collapse = ", "
    )
  )
}

############################################################
# Biotype order
############################################################

preferred_biotype_order_read_barplot <- c(
  "miR",
  "tRF",
  "rRF",
  "snoRNA",
  "piR",
  "ncRNA",
  "cdna"
)

observed_biotypes_read_barplot <- read_gene_status_results %>%
  dplyr::filter(
    !is.na(biotype_priority),
    biotype_priority != "",
    !is.na(feature_id),
    feature_id != ""
  ) %>%
  dplyr::distinct(
    biotype_priority
  ) %>%
  dplyr::pull(
    biotype_priority
  ) %>%
  as.character()

biotype_order_read_barplot <- c(
  intersect(
    preferred_biotype_order_read_barplot,
    observed_biotypes_read_barplot
  ),
  setdiff(
    sort(observed_biotypes_read_barplot),
    preferred_biotype_order_read_barplot
  )
)

if (length(biotype_order_read_barplot) == 0) {
  stop(
    "Aucun biotype n'est disponible pour le barplot read-level."
  )
}

############################################################
# Total number of read-level features per biotype
############################################################

## read_gene_status_results contains one row per feature
## and contrast. Here, keep only one row per
## feature to calculate the total number of distinct features.
total_read_features_by_biotype <- read_gene_status_results %>%
  dplyr::filter(
    !is.na(feature_id),
    feature_id != "",
    !is.na(biotype_priority),
    biotype_priority != ""
  ) %>%
  dplyr::distinct(
    feature_id,
    biotype_priority
  ) %>%
  dplyr::count(
    biotype_priority,
    name = "n_total_features"
  ) %>%
  dplyr::mutate(
    biotype_priority = as.character(
      biotype_priority
    )
  )

print(
  total_read_features_by_biotype
)

############################################################
# Number of up- and downregulated DE read-level features
############################################################

read_de_counts_by_contrast <- read_gene_status_results %>%
  dplyr::filter(
    contrast %in% contrast_order_read_barplot,
    !is.na(feature_id),
    feature_id != "",
    !is.na(biotype_priority),
    biotype_priority != ""
  ) %>%
  dplyr::group_by(
    contrast,
    biotype_priority
  ) %>%
  dplyr::summarise(

    ## All significantly DE read-level features
    n_DE_total = dplyr::n_distinct(
      feature_id[
        read_DE %in% TRUE
      ]
    ),

    ## UP = higher expression in the second group.
    ## Read-level contrasts remain group1 vs group2,
    ## so this corresponds to read_log2FoldChange < 0.
    n_up_DE = dplyr::n_distinct(
      feature_id[
        read_DE %in% TRUE &
          !is.na(read_log2FoldChange) &
          read_log2FoldChange < 0
      ]
    ),

    ## DOWN = higher expression in the first group,
    ## so read_log2FoldChange > 0 in the original DESeq2 table.
    n_down_DE = dplyr::n_distinct(
      feature_id[
        read_DE %in% TRUE &
          !is.na(read_log2FoldChange) &
          read_log2FoldChange > 0
      ]
    ),

    .groups = "drop"
  ) %>%
  dplyr::mutate(
    biotype_priority = as.character(
      biotype_priority
    )
  )

############################################################
# Complete biotype x comparison table
############################################################

read_de_barplot_summary <- tidyr::expand_grid(
  contrast = contrast_order_read_barplot,
  biotype_priority = biotype_order_read_barplot
) %>%
  dplyr::left_join(
    total_read_features_by_biotype,
    by = "biotype_priority"
  ) %>%
  dplyr::left_join(
    read_de_counts_by_contrast,
    by = c(
      "contrast",
      "biotype_priority"
    )
  ) %>%
  dplyr::mutate(

    n_total_features = tidyr::replace_na(
      n_total_features,
      0L
    ),

    n_DE_total = tidyr::replace_na(
      n_DE_total,
      0L
    ),

    n_up_DE = tidyr::replace_na(
      n_up_DE,
      0L
    ),

    n_down_DE = tidyr::replace_na(
      n_down_DE,
      0L
    ),

    n_other_features =
      n_total_features - n_DE_total,

    comparison = factor(
      unname(
        comparison_labels_read_barplot[
          contrast
        ]
      ),
      levels = c(
        "HB",
        "BT",
        "HT"
      )
    ),

    biotype_priority = factor(
      biotype_priority,
      levels = biotype_order_read_barplot
    )
  ) %>%
  dplyr::arrange(
    biotype_priority,
    comparison
  )

############################################################
# Consistency checks
############################################################

if (any(
  read_de_barplot_summary$n_DE_total !=
    read_de_barplot_summary$n_up_DE +
      read_de_barplot_summary$n_down_DE
)) {
  stop(
    "Incohérence : n_DE_total doit être égal à ",
    "n_up_DE + n_down_DE."
  )
}

if (any(
  read_de_barplot_summary$n_DE_total >
    read_de_barplot_summary$n_total_features
)) {
  stop(
    "Incohérence : le nombre de features DE dépasse ",
    "le nombre total de features du biotype."
  )
}

if (any(
  read_de_barplot_summary$n_other_features < 0
)) {
  stop(
    "Incohérence : un nombre négatif de features non-DE ",
    "a été obtenu."
  )
}

print(
  read_de_barplot_summary
)

############################################################
# Export the summary table
############################################################

write.table(
  read_de_barplot_summary,
  file = file.path(
    read_barplot_de_dir,
    paste0(
      "Total_read_features_and_up_down_DE_features_",
      "by_biotype_and_comparison.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Convert to long format
############################################################

read_de_barplot_long <- read_de_barplot_summary %>%
  dplyr::select(
    contrast,
    comparison,
    biotype_priority,
    n_total_features,
    n_other_features,
    n_down_DE,
    n_up_DE
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      n_other_features,
      n_down_DE,
      n_up_DE
    ),
    names_to = "feature_status",
    values_to = "n_features"
  ) %>%
  dplyr::mutate(
    feature_status = factor(
      feature_status,
      levels = c(
        "n_other_features",
        "n_down_DE",
        "n_up_DE"
      ),
      labels = c(
        "Other read-level features",
        "Downregulated DE features",
        "Upregulated DE features"
      )
    )
  )

############################################################
# Readable biotype names
############################################################

biotype_facet_labels_read_barplot <- c(
  miR = "miRNA",
  tRF = "tRF",
  rRF = "rRF",
  snoRNA = "snoRNA",
  piR = "piRNA",
  ncRNA = "ncRNA",
  cdna = "cDNA"
)

############################################################
# Figure
############################################################

output_read_barplot_file <- file.path(
  read_barplot_de_dir,
  "Barplot_total_read_features_and_up_down_DE_by_biotype.pdf"
)

grDevices::pdf(
  output_read_barplot_file,
  width = 13,
  height = 9,
  useDingbats = FALSE
)

p_read_de_barplot_biotype <- ggplot(
  read_de_barplot_long,
  aes(
    x = comparison,
    y = n_features,
    fill = feature_status
  )
) +

  ## The total height corresponds to the total number
  ## of read-level features in the biotype.
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.35
  ) +

  ## Label positioned above the bar.
  ## geom_label preserves readability even for
  ## biotypes containing few features.
  geom_label(
    data = read_de_barplot_summary %>%
      dplyr::filter(
        n_total_features > 0
      ),
    aes(
      x = comparison,
      y = n_total_features,
      label = paste0(
        "Total: ",
        n_total_features,
        "\nUp: ",
        n_up_DE,
        " | Down: ",
        n_down_DE
      )
    ),
    inherit.aes = FALSE,
    vjust = -0.18,
    size = 2.9,
    lineheight = 0.9,
    fill = "white",
    label.size = 0.2
  ) +

  facet_wrap(
    ~ biotype_priority,
    scales = "free_y",
    ncol = 3,
    drop = FALSE,
    labeller = ggplot2::as_labeller(
      biotype_facet_labels_read_barplot,
      default = ggplot2::label_value
    )
  ) +

  scale_fill_manual(
    values = c(
      "Other read-level features" = "grey84",
      "Downregulated DE features" = "#4C6A92",
      "Upregulated DE features" = "#A66A5B"
    ),
    breaks = c(
      "Other read-level features",
      "Upregulated DE features",
      "Downregulated DE features"
    ),
    drop = FALSE
  ) +

  scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.34
      )
    )
  ) +

  labs(
    title = paste0(
      "Total read-level features and up/downregulated ",
      "DE features by biotype"
    ),

    subtitle = paste0(
      "HB = Head vs Body; BT = Body vs Tail; ",
      "HT = Head vs Tail | ",
      "Up: higher in the second group; ",
      "Down: higher in the first group"
    ),

    x = "Comparison",
    y = "Number of distinct read-level features",
    fill = NULL
  ) +

  theme_bw() +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),

    plot.subtitle = element_text(
      hjust = 0.5
    ),

    strip.text = element_text(
      face = "bold"
    ),

    panel.grid.major.x = element_blank(),

    panel.grid.minor = element_blank(),

    legend.position = "top"
  )

print(
  p_read_de_barplot_biotype
)

invisible(
  grDevices::dev.off()
)

message(
  "Barplot read-level up/down terminé."
)

message(
  "Figure générée : ",
  output_read_barplot_file
)


## ----figures-22b-barplot-total-read-features-and-up-down-de-by-biotype, echo=FALSE, message=FALSE, warning=FALSE, results="asis"-----
figs <- c(
  file.path(
    outdir,
    "Barplot_DE_read_features_by_biotype",
    "Barplot_total_read_features_and_up_down_DE_by_biotype.pdf"
  )
)

include_pdf_as_large_png(
  figs
)


## ----25-de-reads-per-trna-gene, eval=TRUE, message=FALSE, warning=FALSE--------------------------------------------------------------
############################################################
# 25. Percentage of DE reads within each tRNA gene
############################################################

message(
  "---- Nombre et pourcentage de reads DE par tRNA ----"
)

trna_read_de_dir <- file.path(
  outdir,
  "DE_reads_within_tRNA_genes"
)

dir.create(
  trna_read_de_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# Check the object produced by block 22
############################################################

if (!exists("per_gene_read_DE_summary")) {
  stop(
    "L'objet per_gene_read_DE_summary est introuvable. ",
    "Le bloc 22 doit être exécuté avant le bloc 25."
  )
}

required_trna_columns <- c(
  "contrast",
  "comparison",
  "gene_id",
  "gene_name_priority",
  "biotype_priority",
  "gene_DE",
  "n_total_reads",
  "n_DE_reads",
  "n_non_DE_reads",
  "percent_DE_reads"
)

missing_trna_columns <- setdiff(
  required_trna_columns,
  colnames(per_gene_read_DE_summary)
)

if (length(missing_trna_columns) > 0) {
  stop(
    "Colonnes manquantes dans per_gene_read_DE_summary : ",
    paste(
      missing_trna_columns,
      collapse = ", "
    )
  )
}

############################################################
# Keep only tRNA / tRF
############################################################

trna_read_de_summary <- per_gene_read_DE_summary %>%
  dplyr::filter(
    as.character(biotype_priority) == "tRF",
    !is.na(gene_name_priority),
    gene_name_priority != ""
  ) %>%
  dplyr::mutate(

    comparison = factor(
      as.character(comparison),
      levels = c(
        "HB",
        "BT",
        "HT"
      )
    ),

    percent_DE_reads = round(
      percent_DE_reads,
      2
    ),

    label_DE_reads = paste0(
      n_DE_reads,
      "/",
      n_total_reads,
      " (",
      sprintf("%.1f", percent_DE_reads),
      "%)"
    )
  )

if (nrow(trna_read_de_summary) == 0) {
  stop(
    "Aucun tRNA/tRF n'a été retrouvé dans ",
    "per_gene_read_DE_summary."
  )
}

############################################################
# Checks
############################################################

if (any(
  trna_read_de_summary$n_DE_reads >
    trna_read_de_summary$n_total_reads
)) {
  stop(
    "Incohérence : un tRNA possède plus de reads DE ",
    "que de reads totaux."
  )
}

if (any(
  trna_read_de_summary$percent_DE_reads < 0 |
    trna_read_de_summary$percent_DE_reads > 100
)) {
  stop(
    "Incohérence : certains pourcentages sont hors ",
    "de l'intervalle 0-100 %."
  )
}

############################################################
# Export the complete table
############################################################

trna_read_de_export <- trna_read_de_summary %>%
  dplyr::select(
    comparison,
    contrast,
    gene_id,
    gene_name_priority,
    n_total_reads,
    n_DE_reads,
    n_non_DE_reads,
    percent_DE_reads,
    gene_DE
  ) %>%
  dplyr::arrange(
    comparison,
    dplyr::desc(percent_DE_reads),
    dplyr::desc(n_DE_reads),
    gene_name_priority
  )

print(
  trna_read_de_export
)

write.table(
  trna_read_de_export,
  file = file.path(
    trna_read_de_dir,
    "DE_reads_percentage_within_each_tRNA.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Overall summary
############################################################

trna_read_de_global_summary <- trna_read_de_summary %>%
  dplyr::group_by(
    comparison
  ) %>%
  dplyr::summarise(

    n_tRNA =
      dplyr::n_distinct(gene_id),

    n_tRNA_with_at_least_one_DE_read =
      dplyr::n_distinct(
        gene_id[n_DE_reads > 0]
      ),

    total_distinct_reads =
      sum(
        n_total_reads,
        na.rm = TRUE
      ),

    total_DE_reads =
      sum(
        n_DE_reads,
        na.rm = TRUE
      ),

    percent_DE_reads_overall =
      dplyr::if_else(
        total_distinct_reads > 0,
        100 *
          total_DE_reads /
          total_distinct_reads,
        0
      ),

    .groups = "drop"
  )

print(
  trna_read_de_global_summary
)

write.table(
  trna_read_de_global_summary,
  file = file.path(
    trna_read_de_dir,
    "Global_DE_reads_summary_tRNA.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# tRNA order in the figure
############################################################

trna_order <- trna_read_de_summary %>%
  dplyr::group_by(
    gene_name_priority
  ) %>%
  dplyr::summarise(

    maximum_percent_DE = max(
      percent_DE_reads,
      na.rm = TRUE
    ),

    maximum_n_DE = max(
      n_DE_reads,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  dplyr::arrange(
    maximum_percent_DE,
    maximum_n_DE,
    gene_name_priority
  ) %>%
  dplyr::pull(
    gene_name_priority
  )

trna_read_de_plot <- trna_read_de_summary %>%
  dplyr::mutate(
    gene_name_priority = factor(
      gene_name_priority,
      levels = trna_order
    )
  )

############################################################
# Automatic figure size
############################################################

n_trna_plot <- dplyr::n_distinct(
  trna_read_de_plot$gene_name_priority
)

trna_pdf_height <- max(
  8,
  0.32 * n_trna_plot + 3
)

############################################################
# Figure
############################################################

output_trna_pdf <- file.path(
  trna_read_de_dir,
  "Percentage_DE_reads_within_each_tRNA.pdf"
)

grDevices::pdf(
  output_trna_pdf,
  width = 18,
  height = trna_pdf_height,
  useDingbats = FALSE
)

p_trna_read_de <- ggplot(
  trna_read_de_plot,
  aes(
    x = percent_DE_reads,
    y = gene_name_priority
  )
) +

  geom_col(
    width = 0.72,
    fill = "#7E93A8",
    color = "white",
    linewidth = 0.3
  ) +

  ## Sufficiently high percentage:
  ## label inside the bar
  geom_text(
    data = trna_read_de_plot %>%
      dplyr::filter(
        percent_DE_reads >= 15
      ),
    aes(
      x = percent_DE_reads - 1,
      label = label_DE_reads
    ),
    hjust = 1,
    color = "white",
    fontface = "bold",
    size = 3
  ) +

  ## Low or zero percentage:
  ## label to the right of the bar
  geom_text(
    data = trna_read_de_plot %>%
      dplyr::filter(
        percent_DE_reads < 15
      ),
    aes(
      x = percent_DE_reads + 1,
      label = label_DE_reads
    ),
    hjust = 0,
    color = "black",
    size = 3
  ) +

  facet_wrap(
    ~ comparison,
    nrow = 1
  ) +

  scale_x_continuous(
    breaks = c(
      0,
      25,
      50,
      75,
      100
    ),
    labels = function(x) {
      paste0(x, "%")
    },
    limits = c(
      0,
      115
    ),
    expand = c(
      0,
      0
    )
  ) +

  labs(
    title = "Percentage of DE reads within each tRNA",
    subtitle = paste0(
      "Labels = number of DE reads / total distinct reads ",
      "assigned to the tRNA"
    ),
    x = "Percentage of reads that are DE",
    y = "tRNA annotation"
  ) +

  theme_bw() +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),

    plot.subtitle = element_text(
      hjust = 0.5
    ),

    strip.text = element_text(
      face = "bold",
      size = 11
    ),

    panel.grid.major.y = element_blank(),

    panel.grid.minor = element_blank(),

    axis.text.y = element_text(
      size = 7
    )
  )

print(
  p_trna_read_de
)

invisible(
  grDevices::dev.off()
)

message(
  "Nombre de tRNA représentés : ",
  n_trna_plot
)

message(
  "Figure générée : ",
  output_trna_pdf
)


## ----figures-25-de-reads-per-trna-gene, echo=FALSE, message=FALSE, warning=FALSE, results="asis"-------------------------------------
figs <- c(
  file.path(
    outdir,
    "DE_reads_within_tRNA_genes",
    "Percentage_DE_reads_within_each_tRNA.pdf"
  )
)

include_pdf_as_large_png(
  figs
)


## ----26-top10-percent-most-abundant-reads, eval=TRUE, message=FALSE, warning=FALSE---------------------------------------------------
############################################################
# 26. Top 10% most abundant read-level features
############################################################

message(
  "---- Top 10 % des read-level features les plus abondantes ----"
)

top10_reads_dir <- file.path(
  outdir,
  "Top10_percent_most_abundant_reads"
)

dir.create(
  top10_reads_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# Check required objects
############################################################

required_objects_top10 <- c(
  "counts_DE",
  "counts_filtered_out",
  "meta",
  "padj_cut",
  "lfc_cut"
)

missing_objects_top10 <- required_objects_top10[
  !vapply(
    required_objects_top10,
    exists,
    logical(1)
  )
]

if (length(missing_objects_top10) > 0) {
  stop(
    "Objets manquants pour le bloc Top 10 % : ",
    paste(
      missing_objects_top10,
      collapse = ", "
    )
  )
}

if (!exists("prepare_spz_metadata", mode = "function")) {
  stop(
    "La fonction prepare_spz_metadata() est introuvable."
  )
}

############################################################
# Read-level feature annotation
############################################################

annotation_top10 <- counts_filtered_out %>%
  dplyr::filter(
    !is.na(feature_id),
    feature_id != ""
  ) %>%
  dplyr::select(
    feature_id,
    sequence,
    gene_name_priority,
    biotype_priority
  ) %>%
  dplyr::distinct(
    feature_id,
    .keep_all = TRUE
  )

############################################################
# Library sizes AFTER filtering
############################################################

filtered_library_sizes <- colSums(
  counts_DE,
  na.rm = TRUE
)

if (any(filtered_library_sizes <= 0)) {
  stop(
    "Au moins un échantillon possède une taille de ",
    "librairie filtrée nulle."
  )
}

############################################################
# Convert read-level counts to CPM
############################################################

counts_DE_cpm <- sweep(
  counts_DE,
  MARGIN = 2,
  STATS = filtered_library_sizes,
  FUN = "/"
) * 1e6

############################################################
# Abundance of each read across the entire filtered dataset
############################################################

read_abundance_all_samples <- tibble::tibble(

  feature_id = rownames(
    counts_DE
  ),

  total_raw_umi_all_samples = rowSums(
    counts_DE,
    na.rm = TRUE
  ),

  mean_raw_umi_all_samples = rowMeans(
    counts_DE,
    na.rm = TRUE
  ),

  mean_cpm_all_samples = rowMeans(
    counts_DE_cpm,
    na.rm = TRUE
  ),

  median_cpm_all_samples = apply(
    counts_DE_cpm,
    1,
    stats::median,
    na.rm = TRUE
  )
) %>%

  dplyr::left_join(
    annotation_top10,
    by = "feature_id"
  )

############################################################
# Exact selection of the top 10% most abundant reads
############################################################

n_filtered_reads <- nrow(
  read_abundance_all_samples
)

n_top10_reads <- max(
  1L,
  ceiling(
    0.10 * n_filtered_reads
  )
)

top10_abundant_reads <- read_abundance_all_samples %>%

  dplyr::arrange(
    dplyr::desc(
      mean_cpm_all_samples
    ),
    dplyr::desc(
      total_raw_umi_all_samples
    ),
    feature_id
  ) %>%

  dplyr::slice_head(
    n = n_top10_reads
  ) %>%

  dplyr::mutate(
    abundance_rank = dplyr::row_number(),

    biotype_label = dplyr::case_when(
      biotype_priority == "miR" ~ "miRNA",
      biotype_priority == "tRF" ~ "tRF",
      biotype_priority == "rRF" ~ "rRF",
      biotype_priority == "snoRNA" ~ "snoRNA",
      biotype_priority == "piR" ~ "piRNA",
      biotype_priority == "ncRNA" ~ "ncRNA",
      biotype_priority == "cdna" ~ "cDNA",

      is.na(biotype_priority) |
        biotype_priority == "" ~
        "Unannotated",

      TRUE ~
        as.character(
          biotype_priority
        )
    )
  )

top10_threshold_cpm <- min(
  top10_abundant_reads$mean_cpm_all_samples,
  na.rm = TRUE
)

message(
  "Nombre total de read-level features après filtration : ",
  n_filtered_reads
)

message(
  "Nombre de features retenues dans le top 10 % : ",
  n_top10_reads
)

message(
  "Seuil minimal de mean CPM dans le top 10 % : ",
  round(
    top10_threshold_cpm,
    4
  )
)

############################################################
# Export the top 10%
############################################################

write.table(
  top10_abundant_reads,
  file = file.path(
    top10_reads_dir,
    "Top10_percent_most_abundant_read_features.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Mean abundance of each read in each major tissue type
############################################################

tissues_top10 <- c(
  "Testis",
  "Epididymis",
  "Spermatozoa"
)

tissues_top10 <- tissues_top10[
  tissues_top10 %in%
    unique(
      as.character(
        meta$`Tissue type`
      )
    )
]

tissue_mean_cpm_top10 <- lapply(
  tissues_top10,
  function(tissue_i) {

    samples_i <- meta %>%
      dplyr::filter(
        as.character(`Tissue type`) ==
          tissue_i
      ) %>%
      dplyr::pull(
        A_sample
      ) %>%
      intersect(
        colnames(counts_DE_cpm)
      )

    if (length(samples_i) == 0) {
      return(NULL)
    }

    tibble::tibble(

      feature_id =
        top10_abundant_reads$feature_id,

      tissue =
        tissue_i,

      mean_cpm_tissue = rowMeans(
        counts_DE_cpm[
          top10_abundant_reads$feature_id,
          samples_i,
          drop = FALSE
        ],
        na.rm = TRUE
      )
    )
  }
) %>%
  dplyr::bind_rows()

if (nrow(tissue_mean_cpm_top10) == 0) {
  stop(
    "Impossible de calculer l'abondance par tissu."
  )
}

############################################################
# Dominant tissue for each read
############################################################

dominant_tissue_top10 <- tissue_mean_cpm_top10 %>%

  dplyr::group_by(
    feature_id
  ) %>%

  dplyr::mutate(
    maximum_tissue_cpm = max(
      mean_cpm_tissue,
      na.rm = TRUE
    ),

    is_maximum =
      abs(
        mean_cpm_tissue -
          maximum_tissue_cpm
      ) < 1e-12
  ) %>%

  dplyr::summarise(

    dominant_tissue = {
      idx_max <- which(
        is_maximum
      )

      if (length(idx_max) == 1) {
        as.character(
          tissue[idx_max]
        )
      } else {
        "Shared"
      }
    },

    .groups = "drop"
  )

top10_abundant_reads <- top10_abundant_reads %>%
  dplyr::left_join(
    dominant_tissue_top10,
    by = "feature_id"
  )

############################################################
# Detailed export of abundance by tissue
############################################################

write.table(
  tissue_mean_cpm_top10,
  file = file.path(
    top10_reads_dir,
    "Top10_percent_reads_mean_CPM_by_tissue.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  top10_abundant_reads,
  file = file.path(
    top10_reads_dir,
    "Top10_percent_reads_with_dominant_tissue.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Pie chart 1: dominant tissue
############################################################

top10_tissue_summary <- top10_abundant_reads %>%

  dplyr::mutate(
    dominant_tissue = dplyr::coalesce(
      dominant_tissue,
      "Unknown"
    )
  ) %>%

  dplyr::count(
    dominant_tissue,
    name = "n_read_features"
  ) %>%

  dplyr::mutate(
    percent =
      100 *
      n_read_features /
      sum(n_read_features)
  ) %>%

  dplyr::arrange(
    dplyr::desc(
      n_read_features
    )
  )

print(
  top10_tissue_summary
)

write.table(
  top10_tissue_summary,
  file = file.path(
    top10_reads_dir,
    "Top10_percent_reads_dominant_tissue_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

tissue_colors_top10 <- c(
  "Testis" = "#607D9D",
  "Epididymis" = "#8FA9C1",
  "Spermatozoa" = "#BCCBD8",
  "Shared" = "#B8B8B8",
  "Unknown" = "#E0E0E0"
)

grDevices::pdf(
  file.path(
    top10_reads_dir,
    "Piechart_top10_percent_reads_by_dominant_tissue.pdf"
  ),
  width = 9,
  height = 7,
  useDingbats = FALSE
)

p_top10_tissue <- ggplot(
  top10_tissue_summary,
  aes(
    x = "",
    y = n_read_features,
    fill = dominant_tissue
  )
) +

  geom_col(
    width = 1,
    color = "white",
    linewidth = 0.6
  ) +

  coord_polar(
    theta = "y"
  ) +

  geom_text(
    aes(
      label = dplyr::if_else(
        percent >= 3,
        paste0(
          round(
            percent,
            1
          ),
          "%"
        ),
        ""
      )
    ),
    position = position_stack(
      vjust = 0.5
    ),
    size = 4
  ) +

  scale_fill_manual(
    values = tissue_colors_top10,
    drop = FALSE
  ) +

  labs(
    title = paste0(
      "Top 10% most abundant reads — ",
      "dominant tissue"
    ),
    subtitle = paste0(
      "One read = one distinct feature | ",
      "dominant tissue = highest mean CPM"
    ),
    fill = "Tissue"
  ) +

  theme_void() +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    legend.position = "right"
  )

print(
  p_top10_tissue
)

invisible(
  grDevices::dev.off()
)

############################################################
# Pie chart 2: biotypes
############################################################

biotype_order_top10 <- c(
  "miRNA",
  "tRF",
  "rRF",
  "snoRNA",
  "piRNA",
  "ncRNA",
  "cDNA",
  "Unannotated"
)

top10_biotype_summary <- top10_abundant_reads %>%

  dplyr::count(
    biotype_label,
    name = "n_read_features"
  ) %>%

  dplyr::mutate(
    percent =
      100 *
      n_read_features /
      sum(n_read_features),

    biotype_label = factor(
      biotype_label,
      levels = c(
        biotype_order_top10,
        setdiff(
          unique(
            as.character(
              biotype_label
            )
          ),
          biotype_order_top10
        )
      )
    )
  ) %>%

  dplyr::arrange(
    biotype_label
  )

print(
  top10_biotype_summary
)

write.table(
  top10_biotype_summary,
  file = file.path(
    top10_reads_dir,
    "Top10_percent_reads_biotype_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

biotype_colors_top10 <- c(
  "miRNA" = "#566D7E",
  "tRF" = "#78909C",
  "rRF" = "#9AAEB8",
  "snoRNA" = "#A9897D",
  "piRNA" = "#8A817C",
  "ncRNA" = "#AAA39D",
  "cDNA" = "#C4B8AE",
  "Unannotated" = "#D9D9D9"
)

grDevices::pdf(
  file.path(
    top10_reads_dir,
    "Piechart_top10_percent_reads_by_biotype.pdf"
  ),
  width = 9,
  height = 7,
  useDingbats = FALSE
)

p_top10_biotype <- ggplot(
  top10_biotype_summary,
  aes(
    x = "",
    y = n_read_features,
    fill = biotype_label
  )
) +

  geom_col(
    width = 1,
    color = "white",
    linewidth = 0.6
  ) +

  coord_polar(
    theta = "y"
  ) +

  geom_text(
    aes(
      label = dplyr::if_else(
        percent >= 3,
        paste0(
          round(
            percent,
            1
          ),
          "%"
        ),
        ""
      )
    ),
    position = position_stack(
      vjust = 0.5
    ),
    size = 3.7
  ) +

  scale_fill_manual(
    values = biotype_colors_top10,
    drop = FALSE
  ) +

  labs(
    title = paste0(
      "Biotype composition of the top 10% ",
      "most abundant reads"
    ),
    subtitle = paste0(
      "Ranking based on mean CPM across all ",
      "filtered samples"
    ),
    fill = "Biotype"
  ) +

  theme_void() +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    legend.position = "right"
  )

print(
  p_top10_biotype
)

invisible(
  grDevices::dev.off()
)

############################################################
# Differential analysis of SPZ Head vs Tail
# restricted to the top 10%
############################################################

message(
  "---- DESeq2 Top 10 % : SPZ Head vs Tail ----"
)

meta_top10_spz_ht <- prepare_spz_metadata(
  regions = c(
    "Head",
    "Tail"
  )
)

meta_top10_spz_ht$Origin_in_Epididymis <- relevel(
  meta_top10_spz_ht$Origin_in_Epididymis,
  ref = "Head"
)

top10_spz_samples <- intersect(
  rownames(
    meta_top10_spz_ht
  ),
  colnames(
    counts_DE
  )
)

meta_top10_spz_ht <- meta_top10_spz_ht[
  top10_spz_samples,
  ,
  drop = FALSE
]

counts_top10_spz_ht <- counts_DE[
  top10_abundant_reads$feature_id,
  top10_spz_samples,
  drop = FALSE
]

############################################################
# Remove only features absent from all SPZ samples
############################################################

features_present_top10_spz <- rowSums(
  counts_top10_spz_ht,
  na.rm = TRUE
) > 0

counts_top10_spz_ht <- counts_top10_spz_ht[
  features_present_top10_spz,
  ,
  drop = FALSE
]

message(
  "Features du Top 10 % présentes dans les SPZ Head/Tail : ",
  nrow(
    counts_top10_spz_ht
  ),
  " / ",
  n_top10_reads
)

if (nrow(counts_top10_spz_ht) < 2) {
  stop(
    "Moins de deux features du Top 10 % sont présentes ",
    "dans les spermatozoïdes Head/Tail."
  )
}

############################################################
# DESeq2
############################################################

dds_top10_spz_ht <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(
    counts_top10_spz_ht
  ),
  colData = meta_top10_spz_ht,
  design = ~ Origin_in_Epididymis
)

dds_top10_spz_ht <- DESeq2::DESeq(
  dds_top10_spz_ht
)

res_top10_spz_ht <- DESeq2::results(
  dds_top10_spz_ht,
  contrast = c(
    "Origin_in_Epididymis",
    "Head",
    "Tail"
  )
)

############################################################
# Annotate results
############################################################

res_top10_spz_ht_df <- as.data.frame(
  res_top10_spz_ht
) %>%

  tibble::rownames_to_column(
    "feature_id"
  ) %>%

  dplyr::left_join(
    top10_abundant_reads %>%
      dplyr::select(
        feature_id,
        abundance_rank,
        mean_cpm_all_samples,
        total_raw_umi_all_samples,
        sequence,
        gene_name_priority,
        biotype_priority,
        biotype_label,
        dominant_tissue
      ),
    by = "feature_id"
  ) %>%

  dplyr::mutate(

    neglog10padj = dplyr::if_else(
      is.na(padj),
      NA_real_,
      -log10(
        pmax(
          padj,
          .Machine$double.xmin
        )
      )
    ),

    ## The DESeq2 contrast remains Head vs Tail.
    ## Create a variable for display only so that
    ## Tail appears on the right side of the volcano plot without changing the original log2FC.
    plot_log2FoldChange = -log2FoldChange,

    signif =
      !is.na(padj) &
      padj < padj_cut &
      abs(log2FoldChange) >= lfc_cut,

    direction = dplyr::case_when(

      signif &
        log2FoldChange > 0 ~
        "Higher in Head",

      signif &
        log2FoldChange < 0 ~
        "Higher in Tail",

      TRUE ~
        "Not DE"
    )
  )

############################################################
# DE summary
############################################################

n_top10_tested <- nrow(
  res_top10_spz_ht_df
)

n_top10_with_padj <- sum(
  !is.na(
    res_top10_spz_ht_df$padj
  )
)

n_top10_DE <- sum(
  res_top10_spz_ht_df$signif,
  na.rm = TRUE
)

n_top10_higher_head <- sum(
  res_top10_spz_ht_df$direction ==
    "Higher in Head",
  na.rm = TRUE
)

n_top10_higher_tail <- sum(
  res_top10_spz_ht_df$direction ==
    "Higher in Tail",
  na.rm = TRUE
)

top10_DE_summary <- tibble::tibble(
  n_filtered_read_features =
    n_filtered_reads,

  n_top10_features =
    n_top10_reads,

  n_top10_features_present_in_SPZ =
    n_top10_tested,

  n_features_with_padj =
    n_top10_with_padj,

  n_DE_features =
    n_top10_DE,

  n_higher_in_Head =
    n_top10_higher_head,

  n_higher_in_Tail =
    n_top10_higher_tail
)

print(
  top10_DE_summary
)

write.table(
  top10_DE_summary,
  file = file.path(
    top10_reads_dir,
    "Top10_percent_reads_SPZ_Head_vs_Tail_DE_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  res_top10_spz_ht_df,
  file = file.path(
    top10_reads_dir,
    "DESeq2_Top10_percent_reads_SPZ_Head_vs_Tail.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# Features to label
############################################################

top_labels_top10_volcano <- res_top10_spz_ht_df %>%

  dplyr::filter(
    signif,
    !is.na(padj)
  ) %>%

  dplyr::arrange(
    padj
  ) %>%

  dplyr::slice_head(
    n = 15
  ) %>%

  dplyr::mutate(

    label_plot = dplyr::case_when(

      !is.na(gene_name_priority) &
        gene_name_priority != "" ~

        paste0(
          gene_name_priority,
          " [#",
          abundance_rank,
          "]"
        ),

      !is.na(sequence) &
        sequence != "" ~

        paste0(
          substr(
            sequence,
            1,
            15
          ),
          "... [#",
          abundance_rank,
          "]"
        ),

      TRUE ~

        paste0(
          "feature #",
          abundance_rank
        )
    )
  )

############################################################
# Volcano plot
############################################################

volcano_colors_top10 <- c(
  "Not DE" = "#C7C7C7",
  "Higher in Head" = "#617C93",
  "Higher in Tail" = "#A77768"
)

grDevices::pdf(
  file.path(
    top10_reads_dir,
    "Volcano_Top10_percent_reads_SPZ_Head_vs_Tail.pdf"
  ),
  width = 11,
  height = 9,
  useDingbats = FALSE
)

p_top10_volcano <- ggplot(
  res_top10_spz_ht_df,
  aes(
    x = plot_log2FoldChange,
    y = neglog10padj,
    color = direction
  )
) +

  geom_point(
    alpha = 0.65,
    size = 1.5
  ) +

  geom_vline(
    xintercept = c(
      -lfc_cut,
      lfc_cut
    ),
    linetype = "dashed",
    linewidth = 0.5
  ) +

  geom_hline(
    yintercept = -log10(
      padj_cut
    ),
    linetype = "dashed",
    linewidth = 0.5
  ) +

  ggrepel::geom_text_repel(
    data = top_labels_top10_volcano,
    aes(
      label = label_plot
    ),
    color = "black",
    size = 3,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.25,
    min.segment.length = 0
  ) +

  annotate(
    "label",
    x = Inf,
    y = Inf,

    label = paste0(
      "Top 10% selected: ",
      n_top10_reads,
      "\nTested in SPZ: ",
      n_top10_tested,
      "\nDE: ",
      n_top10_DE,
      "\nHigher in Head: ",
      n_top10_higher_head,
      "\nHigher in Tail: ",
      n_top10_higher_tail
    ),

    hjust = 1.02,
    vjust = 1.02,

    size = 3.4,
    label.size = 0.25,
    fill = "white",
    alpha = 0.9
  ) +

  scale_color_manual(
    values = volcano_colors_top10,
    breaks = c(
      "Not DE",
      "Higher in Head",
      "Higher in Tail"
    )
  ) +

  labs(
    title = paste0(
      "Top 10% most abundant reads — ",
      "SPZ Head vs Tail"
    ),

    subtitle = paste0(
      "Read-level DESeq2 | Top 10% selected by mean CPM | ",
      "padj < ",
      padj_cut,
      " and |log2FC| ≥ ",
      lfc_cut
    ),

    x = "Displayed log2(FoldChange): Tail / Head",

    y = "-log10(padj)",

    color = NULL
  ) +

  theme_bw() +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),

    plot.subtitle = element_text(
      hjust = 0.5
    ),

    panel.grid.minor = element_blank(),

    legend.position = "top"
  )

print(
  p_top10_volcano
)

invisible(
  grDevices::dev.off()
)

message(
  "Analyse du Top 10 % terminée."
)

message(
  "Résultats dans : ",
  top10_reads_dir
)


## ----figures-26-top10-percent-most-abundant-reads, echo=FALSE, message=FALSE, warning=FALSE, results="asis"--------------------------
figs <- c(

  file.path(
    outdir,
    "Top10_percent_most_abundant_reads",
    "Piechart_top10_percent_reads_by_dominant_tissue.pdf"
  ),

  file.path(
    outdir,
    "Top10_percent_most_abundant_reads",
    "Piechart_top10_percent_reads_by_biotype.pdf"
  ),

  file.path(
    outdir,
    "Top10_percent_most_abundant_reads",
    "Volcano_Top10_percent_reads_SPZ_Head_vs_Tail.pdf"
  )
)

include_pdf_as_large_png(
  figs
)


## ----check-all-zero-lines-before-filtering, eval=TRUE---------------------------------------------------------------------------------
############################################################
# Check rows containing only zeros before filtering
############################################################

message("---- Vérification lignes avec uniquement des 0 avant filtration ----")

## counts_only = raw count matrix, before counts_filt
zero_rows_before_filter <- rowSums(counts_only, na.rm = TRUE) == 0

n_zero_rows_before_filter <- sum(zero_rows_before_filter)

message("Nombre de lignes avec uniquement des 0 avant filtration : ",
        n_zero_rows_before_filter)

message("Pourcentage : ",
        round(100 * n_zero_rows_before_filter / nrow(counts_only), 3),
        "%")

## View the affected rows in the complete table with annotations
zero_rows_table <- counts[zero_rows_before_filter, , drop = FALSE]

print(head(zero_rows_table, 20))

write.table(
  zero_rows_table,
  file = file.path(outdir, "rows_all_zero_before_filter.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
