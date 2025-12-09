#!/bin/bash
#BSUB -J nx120
#BSUB -o logs/heart_mpra_%J.out
#BSUB -e logs/heart_mpra_%J.err
#BSUB -M 60GB
#BSUB -n 20
#BSUB -R "rusage[mem=60GB] span[hosts=1]"
#BSUB -W 24:00

# Load required modules
module load apptainer
conda activate snakemake 

# Set working directory and create logs directory
WORKDIR="/home/fengyq/data/Ning_MPRA/RNA_DNA/count_ratio/NX120"
cd "${WORKDIR}"
mkdir -p logs

# Configuration files
CONFIG="${WORKDIR}/Di_mpra_experiment.yaml"
SNAKEFILE="/home/fengyq/data/Ning_MPRA/MPRAsnakeflow/workflow/Snakefile"


# Run the workflow
snakemake \
    --configfile "${CONFIG}" \
    --snakefile "${SNAKEFILE}" \
    --sdm conda \
    --notemp \
    --cores 20 \
    --rerun-incomplete

