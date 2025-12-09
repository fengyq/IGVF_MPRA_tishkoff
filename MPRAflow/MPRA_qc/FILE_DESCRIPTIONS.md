# MPRA_qc Directory File Descriptions

## Quality Control Reports

### `MPRAsnakeflow_assignment.QC_report.html`
**Purpose**: Comprehensive HTML quality control report for the assignment workflow

## Quality Metrics Visualizations

### `NX120_nBC1_dna_normalized_cor_heatmap.pdf`
**Purpose**: Heatmap of DNA barcode correlations between replicates (normalized)
- Shows correlation matrix for DNA input counts across samples
- Minimum 1 barcode threshold (nBC1) applied
- Color-coded correlation values (high = red, low = blue)
- Used to assess technical reproducibility of DNA library preparation
- Expected high correlations (>0.9) between biological replicates

### `NX120_nBC1_rna_normalized_cor_heatmap.pdf`
**Purpose**: Heatmap of RNA barcode correlations between replicates (normalized)
- Shows correlation matrix for RNA output counts across samples
- Minimum 1 barcode threshold (nBC1) applied
- Reveals biological variability in regulatory activity
- Used to assess consistency of regulatory measurements
- Expected moderate-high correlations (>0.8) for consistent regulation

### `NX120_nBC20_log2FC_violin_by_label.pdf`
**Purpose**: Violin plot of log2 fold changes across all oligos (20 barcode minimum)
- Shows distribution of regulatory effects for each oligos
- Minimum 20 barcode threshold (nBC20) for statistical confidence

### `NX120.Hadza_neg_pos_n_bc_hist.pdf`
**Purpose**: Histogram of barcode counts for MPRA oligos