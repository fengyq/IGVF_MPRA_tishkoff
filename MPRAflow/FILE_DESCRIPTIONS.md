# MPRAflow File Descriptions

## Configuration Files

### `Di_mpra_assignment.yaml`
**Purpose**: Configuration file for the barcode assignment workflow
- Defines barcode structure (15 bp length)
- Specifies adapter sequences for read trimming
- Sets BWA mapping parameters (min_mapping_quality: 1)
- Configures NGmerge parameters for read merging
- Contains input file paths for sequencing data
- Defines filtering thresholds (min_support: 3, fraction: 0.75)

### `Di_mpra_experiment.yaml`
**Purpose**: Configuration file for the experiment-level analysis workflow
- Defines experiment parameters (bc_length: 15, umi_length: 10)
- Contains data folder and experiment file paths
- Specifies assignment file location from previous stage
- Includes variant mapping file path
- Sets count thresholds (bc_threshold: 20, min_counts: 1-2)
- Configures barcode count requirements for replicates (40, 60)

### `heart_mpra_experiment_NX120.csv`
**Purpose**: Sample sheet for MPRA experiment
- Lists sample metadata and file paths
- Contains sample names, conditions, and replicate information
- Used by experiment workflow to locate input FASTQ files

## Reference Data Files

### `Di_mpra_igvfupload.heart.ref_oligos.fa`
**Purpose**: Reference oligo sequences for heart MPRA library
- Contains REF and ALT allele sequences in FASTA format
- Headers include ID, REF, and ALT information for variant mapping

### `Di_mpra_igvfupload.heart.ref_oligos.label.tsv`
**Purpose**: Label file containing metadata for each oligo

## Workflow Scripts

### `Di_mpra_assignment_lsf.sh`
**Purpose**: LSF submission script for assignment workflow on HPC

### `Di_mpra_experiment_lsf.sh`
**Purpose**: LSF submission script for experiment workflow on HPC

## Workflow Rules

### `mapping_bwa_revised.smk`
**Purpose**: Modified Snakemake rules for BWA-based barcode assignment
- Contains rules for reference preparation and indexing
- Implements read mapping with BWA MEM
- Includes enhanced barcode extraction with strict filters
- Rules: assignment_mapping_bwa_ref, assignment_mapping_bwa, assignment_mapping_bwa_getBCs
- Applies zero-mismatch (NM=0) and pure-match CIGAR filters
- Handles BAM merging and indexing for final output

## Documentation

### `0_MPRAsnakeflow_diagram.md`
**Purpose**: Detailed ASCII diagram of the MPRA analysis workflow

## Visualization Files

### `mpra.workflow.pdf`
**Purpose**: Complete workflow DAG (Directed Acyclic Graph) visualization

### `mpra.rulegraph.pdf`
**Purpose**: Rule graph showing workflow logic
