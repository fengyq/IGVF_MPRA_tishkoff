#!/usr/bin/env Rscript

# ===============================================================
# Author: fengyq
# Date: 2025-01-16
# Description: MPRA Analysis - Test Elements vs Negative Controls
#              Performs differential activity analysis using mpralm
#              and generates publication-ready violin plots with
#              statistical comparisons using ggpubr
# 
# Input files:
#   - rawcounts/reporter_experiment.oligo.NX120.fromFile.default.all_rawcount.tsv
#   - heart_CRE_uniq.label_test_control.tsv
# 
# Output:
#   - Volcano plot: NX120_mpralm_volcano.pdf
#   - Violin plots: NX120_logFC_violin.pdf, NX120_logFC_violin_enhanced.pdf
#   - Results table: NX120_test_vs_negctrl_mpralm.tsv
# 
# Usage: Rscript mpra_test_vs_negctrl.R
# ===============================================================

# ===============================================================
# SECTION 1: SETUP AND CONFIGURATION
# ===============================================================

# Load required libraries
suppressPackageStartupMessages({
  library(mpra)        # MPRA analysis framework
  library(tidyverse)   # Data manipulation and visualization
  library(EnhancedVolcano)  # Enhanced volcano plots
  library(ggpubr)      # Publication-ready plots
  library(rstatix)     # Statistical testing
})

# Configuration: Experiment parameters and file paths
experiment_id <- "NX120"  # Change this for different experiments (e.g., "NX118", "NX119", "NX120")

# Input files (automatically constructed based on experiment_id)
counts_file <- paste0("rawcounts/reporter_experiment.oligo.", experiment_id, ".fromFile.default.all_rawcount.tsv")
group_file  <- "heart_CRE_uniq.label_test_control.tsv"

# Output directory and file naming convention
plot_dir    <- "CRE_activity"  # Directory for all CRE activity analysis outputs
output_prefix <- paste0(experiment_id, "_cre_activity")  # Prefix for all output files

# Generated file names will be:
# - Results: {experiment_id}_cre_activity_test_vs_negctrl_mpralm.tsv
# - Volcano: {experiment_id}_cre_activity_mpralm_volcano.pdf  
# - Violin:  {experiment_id}_cre_activity_logFC_violin.pdf

# Set options for cleaner output
options(na.print = "NA")

# Print configuration for verification
message("[config] experiment_id: ", experiment_id)
message("[config] counts_file:   ", counts_file)
message("[config] group_file:    ", group_file)
message("[config] output_prefix: ", output_prefix)
message("[config] plot_dir:      ", plot_dir)

# ===============================================================
# SECTION 2: DATA LOADING AND PREPROCESSING
# ===============================================================

# Load raw count data (DNA and RNA counts per oligo per replicate)
counts <- readr::read_tsv(
  counts_file,
  col_types = readr::cols(
    replicate   = readr::col_integer(),
    oligo_name  = readr::col_character(),
    dna_counts  = readr::col_double(),
    rna_counts  = readr::col_double()
  )
)

# Validate that count data was loaded successfully
if (nrow(counts) == 0) stop("Counts table is empty: ", counts_file)

# Load group assignments (test elements vs negative controls)
groups <- readr::read_tsv(
  group_file,
  col_names = c("oligo_name", "group"),
  col_types = readr::cols(
    oligo_name = readr::col_character(),
    group      = readr::col_factor(levels = c("NegCtrl", "test"))
  )
) %>%
  dplyr::distinct(oligo_name, .keep_all = TRUE)

# Merge count data with group assignments
counts_labeled <- counts %>%
  dplyr::inner_join(groups, by = "oligo_name")

# Validate that we have overlapping data
if (nrow(counts_labeled) == 0) stop("No overlapping oligos between counts and group labels.")

# Identify available replicates
rep_levels <- sort(unique(counts_labeled$replicate))
message("[config] replicates detected: ", paste(rep_levels, collapse = ", "))

# Complete the data matrix (fill missing combinations with NA)
counts_complete <- counts_labeled %>%
  tidyr::complete(oligo_name, replicate = rep_levels) %>%
  dplyr::arrange(oligo_name, replicate)

dna_wide <- counts_complete %>%
  dplyr::select(oligo_name, replicate, dna_counts) %>%
  tidyr::pivot_wider(names_from = replicate, values_from = dna_counts, names_prefix = "rep")

rna_wide <- counts_complete %>%
  dplyr::select(oligo_name, replicate, rna_counts) %>%
  tidyr::pivot_wider(names_from = replicate, values_from = rna_counts, names_prefix = "rep")

# ===============================================================
# SECTION 3: MPRA ANALYSIS USING MPRALM
# ===============================================================

# Ensure DNA and RNA matrices have matching row order
stopifnot(identical(dna_wide$oligo_name, rna_wide$oligo_name))

# Convert to matrices for MPRA analysis (rows = oligos, columns = replicates)
dna_mat <- dna_wide %>% tibble::column_to_rownames("oligo_name") %>% as.matrix()
rna_mat <- rna_wide %>% tibble::column_to_rownames("oligo_name") %>% as.matrix()

# Print matrix dimensions for verification
message("[matrix] DNA dims: ", paste(dim(dna_mat), collapse = " x "))
message("[matrix] RNA dims: ", paste(dim(rna_mat), collapse = " x "))

# Create MPRASet object for analysis
mpra_obj <- MPRASet( DNA = dna_mat,
  RNA = rna_mat,
  eid = rownames(dna_mat),
  eseq = NULL,
  barcode = NULL
)

# Create design matrix (intercept-only model for activity estimation)
design <- matrix(1, nrow = ncol(dna_mat), ncol = 1,
                 dimnames = list(colnames(dna_mat), "intcpt"))

# Fit MPRA linear model to estimate element activities
fit <- mpralm(
  object = mpra_obj,
  design = design,
  aggregate = "mean",
  normalize = TRUE,
  model_type = "indep_groups",
  plot = TRUE
)

top <- limma::topTable(fit, coef = 1, number = Inf) %>%
  tibble::rownames_to_column(var = "oligo_name") %>%
  dplyr::left_join(groups, by = "oligo_name")

# Validate group assignment
if (any(is.na(top$group))) warning("Some elements missing group labels after join.")

# Extract logFC values from negative controls to set activity threshold
neg_logfc <- top %>%
  dplyr::filter(group == "NegCtrl") %>%
  dplyr::pull(logFC)

# Ensure we have negative controls for threshold calculation
if (length(neg_logfc) == 0) stop("No negative controls found to set threshold.")

# Set activity threshold as 95th percentile of negative control distribution
neg_threshold <- stats::quantile(neg_logfc, probs = 0.95, na.rm = TRUE)

# Call active elements based on threshold and significance
activity_calls <- top %>%
  dplyr::mutate(
    active = dplyr::case_when(
      group == "test" & logFC > neg_threshold & adj.P.Val < 0.05 ~ TRUE,
      TRUE ~ FALSE
    ),
    activity_vs_neg = logFC - neg_threshold,
    neg_threshold = neg_threshold
  )

summary_tbl <- activity_calls %>%
  dplyr::arrange(dplyr::desc(active), adj.P.Val) %>%
  dplyr::select(oligo_name, group, logFC, adj.P.Val, active, activity_vs_neg) %>%
  head(20)

cat("\n[threshold] NegCtrl 95th percentile logFC:", neg_threshold, "\n")
cat("[summary] Active test elements:", sum(activity_calls$active, na.rm = TRUE), "\n\n")

summary_tbl_df <- as.data.frame(summary_tbl)
print(summary_tbl_df, row.names = FALSE)

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

result_path <- paste0(output_prefix, "_test_vs_negctrl_mpralm.tsv")
readr::write_tsv(activity_calls, result_path)
cat("\n[output] Full results written to:", result_path, "\n")

# ===============================================================
# SECTION 5: PUBLICATION-READY VISUALIZATIONS
# ===============================================================

# ---- Enhanced Volcano Plot ----
# Prepare data for volcano plot (requires rownames)
top_for_volcano <- activity_calls %>% tibble::column_to_rownames("oligo_name")

# Select significant elements with large effect sizes for labeling
selectLab <- rownames(top_for_volcano)[activity_calls$active & abs(activity_calls$logFC) > 1]

# Create enhanced volcano plot
p2 <- EnhancedVolcano(
  top_for_volcano,
  lab = rownames(top_for_volcano),
  selectLab = selectLab,
  x = "logFC",
  y = "adj.P.Val",
  borderWidth = 0.25,
  labSize = 2,
  title = NULL,
  titleLabSize = 16,
  subtitle = NULL,
  axisLabSize = 16,
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 1,
  colAlpha = 0.8,
  shape = 19,
  legendLabSize = 15,
  captionLabSize = 15,
  cutoffLineType = "longdash",
  cutoffLineWidth = 0.25,
  drawConnectors = TRUE
)

print(p2)

volcano_file <- file.path(plot_dir, paste0(output_prefix, "_mpralm_volcano.pdf"))

ggplot2::ggsave( filename = volcano_file, plot = p2, width = 6, height = 5,units = "in", bg = "white")

cat("[plot] Volcano saved to:", volcano_file, "\n")

# ---- Publication-ready Violin Plot of logFC by group ----

# Perform statistical test for comparison
stat_test <- activity_calls %>%
  rstatix::wilcox_test(logFC ~ group) %>%
  rstatix::add_significance() %>%
  rstatix::add_xy_position(x = "group")

# Create publication-ready violin plot
violin_plot <- ggviolin(
  activity_calls, 
  x = "group", 
  y = "logFC",
  fill = "group",
  palette = c("NegCtrl" = "#E7E7E7", "test" = "#2E86AB"),
  alpha = 0.7,
  trim = FALSE,
  draw_quantiles = c(0.25, 0.5, 0.75),
  add = "boxplot",
  add.params = list(fill = "white", alpha = 0.8, width = 0.1)
) +
  geom_hline(yintercept = neg_threshold, linetype = "dashed", 
             colour = "grey30", linewidth = 0.8) +
  stat_pvalue_manual(stat_test, 
                     label = "p.signif", 
                     tip.length = 0.01,
                     bracket.size = 0.6,
                     textsize = 5) +
  labs(
    title = "MPRA Activity Distribution",
    subtitle = "Comparison between Test Elements and Negative Controls",
    x = "Element Type",
    y = "log2(RNA/DNA) Activity Score",
    caption = paste0("Dashed line: negative threshold (", round(neg_threshold, 2), ")")
  ) +
  theme_pubr(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    plot.caption = element_text(size = 9, color = "grey50"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  scale_x_discrete(labels = c("NegCtrl" = "Negative\nControl", "test" = "Test\nElements")) +
  coord_cartesian(clip = "off")

violin_file <- file.path(plot_dir, paste0(output_prefix, "_logFC_violin.pdf"))

ggplot2::ggsave(
  filename = violin_file,
  plot = violin_plot,
  width = 5.5,
  height = 5,
  units = "in",
  bg = "white"
)

cat("[plot] Violin saved to:", violin_file, "\n")
