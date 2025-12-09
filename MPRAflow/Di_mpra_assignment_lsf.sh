#!/bin/bash
#BSUB -J mpra1
#BSUB -o logs/heart_mpra_%J.out
#BSUB -e logs/heart_mpra_%J.err
#BSUB -M 20GB
#BSUB -n 20
#BSUB -R "rusage[mem=20GB] span[hosts=1]"
#BSUB -W 24:00

# Load required modules
module load apptainer
conda activate snakemake

# Set working directory and create logs directory
WORKDIR="/home/fengyq/data/Ning_MPRA/association"
cd "${WORKDIR}"
mkdir -p logs

# Configuration files
CONFIG="${WORKDIR}/heart_mpra_assignment.yaml"
SNAKEFILE="/home/fengyq/data/Ning_MPRA/MPRAsnakeflow/workflow/Snakefile"

# Generate workflow diagram
snakemake --configfile "${CONFIG}" \
    --snakefile "${SNAKEFILE}" \
    --dag | dot -Tpdf -o mpra.workflow.pdf

snakemake --configfile "${CONFIG}" \
    --snakefile "${SNAKEFILE}" \
    --rulegraph | dot -Tpdf -o mpra.rulegraph.pdf

# Generate workflow pseudocode
snakemake --dry-run -p \
    --configfile "${CONFIG}" \
    --snakefile "${SNAKEFILE}" \
    > workflow.pseudocode.txt 2>&1

# Run the workflow
snakemake \
    --configfile "${CONFIG}" \
    --snakefile "${SNAKEFILE}" \
    --sdm conda \
    --notemp \
    --cores 20 \
    --rerun-incomplete

# Generate workflow summary
snakemake --summary \
    --configfile "${CONFIG}" \
    --snakefile "${SNAKEFILE}" \
    > workflow.summary.tsv


# then delete tmp
# snakemake --delete-temp-output
