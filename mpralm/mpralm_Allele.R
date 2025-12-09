## ----------------- MPRA allele comparison (NX120) -----------------
# Compare DNA/RNA ratios between Ref and Alt alleles across 4 replicates
# Inputs (in working directory):
#   - reporter_experiment.oligo.NX120.fromFile.default.all_rawcount.tsv
#   - heart_variants_map.tsv (columns: ID, REF, ALT; values are oligo names)

suppressPackageStartupMessages({
  library(mpra)
  library(tidyverse)
  library(ggrepel)
  library(ggpubr)
  library(ggExtra)
  library(reshape2)
  library(EnhancedVolcano)
  library(scales)
})

# ---------- File paths ----------
counts_file <- "reporter_experiment.oligo.NX120.fromFile.default.all_rawcount.tsv"
map_file    <- "heart_variants_map.tsv"

# ---------- Read data ----------
counts <- readr::read_tsv(
  counts_file,
  col_types = readr::cols(
    replicate   = readr::col_integer(),
    oligo_name  = readr::col_character(),
    dna_counts  = readr::col_double(),
    rna_counts  = readr::col_double()
  )
)

variant_map <- readr::read_tsv(
  map_file,
  col_types = readr::cols(
    ID  = readr::col_character(),
    REF = readr::col_character(),
    ALT = readr::col_character()
  )
)

# ---------- Label each oligo as ref/alt and attach variant ID ----------
labels_long <- variant_map %>%
  tidyr::pivot_longer(c(REF, ALT), names_to = "allele", values_to = "oligo_name") %>%
  dplyr::mutate(allele = tolower(allele))

counts_labeled <- counts %>%
  dplyr::inner_join(labels_long, by = "oligo_name") %>%
  dplyr::select(replicate, ID, allele, oligo_name, dna_counts, rna_counts)

# ---------- Choose variants with ≥3 replicates having both alleles ----------
rep_levels <- sort(unique(counts_labeled$replicate))

id_rep_alleles <- counts_labeled %>%
  dplyr::group_by(ID, replicate) %>%
  dplyr::summarise(n_alleles = dplyr::n_distinct(allele), .groups = "drop")

complete_reps <- id_rep_alleles %>% dplyr::filter(n_alleles == 2L)

valid_ids <- complete_reps %>%
  dplyr::count(ID, name = "n_complete_reps") %>%
  dplyr::filter(n_complete_reps >= 3L) %>%
  dplyr::pull(ID)

# Keep only replicates where both alleles observed; later complete to fill missing replicate as NA
counts_valid <- counts_labeled %>%
  dplyr::semi_join(complete_reps, by = c("ID", "replicate")) %>%
  dplyr::filter(ID %in% valid_ids)

# ---------- Build DNA/RNA matrices for mpralm (columns: ref_1..4, alt_1..4) ----------
counts_for_mpralm <- counts_valid %>%
  dplyr::group_by(ID) %>%
  tidyr::complete(replicate = rep_levels, allele = c("ref", "alt")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(colname = paste0(allele, "_", replicate))

dna_wide <- counts_for_mpralm %>%
  dplyr::select(ID, colname, dna_counts) %>%
  tidyr::pivot_wider(names_from = colname, values_from = dna_counts)

rna_wide <- counts_for_mpralm %>%
  dplyr::select(ID, colname, rna_counts) %>%
  tidyr::pivot_wider(names_from = colname, values_from = rna_counts)

# Ensure consistent column order
col_name <- c(paste0("ref_", rep_levels), paste0("alt_", rep_levels))

Di_DNA <- dna_wide %>% tibble::column_to_rownames("ID") %>% as.data.frame()
Di_RNA <- rna_wide %>% tibble::column_to_rownames("ID") %>% as.data.frame()

# Add any missing columns and order
missing_cols_dna <- setdiff(col_name, colnames(Di_DNA))
missing_cols_rna <- setdiff(col_name, colnames(Di_RNA))
if (length(missing_cols_dna) > 0) Di_DNA[, missing_cols_dna] <- NA_real_
if (length(missing_cols_rna) > 0) Di_RNA[, missing_cols_rna] <- NA_real_
Di_DNA <- Di_DNA[, col_name]
Di_RNA <- Di_RNA[, col_name]

# Column sums for quick sanity check
print(colSums(Di_DNA, na.rm = TRUE))
print(colSums(Di_RNA, na.rm = TRUE))

## -------- construct a MPRASet and run mpralm (paired 4 replicates) --------
mpra_Di <- MPRASet(DNA = as.matrix(Di_DNA), RNA = as.matrix(Di_RNA), eid = rownames(Di_DNA), eseq = NULL, barcode = NULL)

# the design of the MPRASet: 4 ref, 4 alt
design <- data.frame(intcpt = 1, alleleB = grepl("alt", colnames(mpra_Di)))
# block is a vector that pairs columns across replicates
block_vector <- rep(seq_along(rep_levels), 2)

# Linear models for differential analysis of MPRA data
mpralm_Di_fit <- mpralm(object = mpra_Di, design = design, aggregate = "none", normalize = TRUE, block = block_vector, model_type = "corr_groups", plot = TRUE)
toptab_Di <- topTable(mpralm_Di_fit, coef = 2, number = Inf)
head(toptab_Di)
toptab_Di_ranked <- toptab_Di %>% tibble::rownames_to_column(var = "ID") %>% dplyr::arrange(adj.P.Val)

# Write the mpralm results
write.table(toptab_Di_ranked, file = "NX120_Ref_vs_Alt_mpralm_sigdiff.tsv", row.names = FALSE, sep = "\t", quote = FALSE)

# ---------------------------------------- Volcano plot ----------------------------------------
selectLab2 <- toptab_Di %>% dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 0.25) %>% rownames()
print(summary(-log10(toptab_Di$adj.P.Val)))
print(summary(toptab_Di$logFC))
p2 <- EnhancedVolcano( toptab_Di, lab = rownames(toptab_Di),selectLab = selectLab2, x = "logFC", y = "adj.P.Val", borderWidth = 0.25, labSize = 2, title = NULL,  titleLabSize = 16, subtitle = NULL, axisLabSize = 16, pCutoff = 0.05,  FCcutoff = 1,  pointSize = 2,  colAlpha = 0.8, shape = 19, legendLabSize = 15, captionLabSize = 15, cutoffLineType = 'longdash', cutoffLineWidth = 0.25, drawConnectors = TRUE)
print(p2)
ggsave("NX120_Ref_vs_Alt_mpralm_sigdiff.Volcano.pdf", width = 5, height = 6, units = "in")
