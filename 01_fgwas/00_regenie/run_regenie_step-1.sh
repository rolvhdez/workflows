#!/bin/bash

# Author: Roberto Olvera Hernandez
# Date: 2026-03-02
#
# --------------------------------
# REGENIE STEP 01
# --------------------------------
# Calculate null model (whole-genome regression)
#
# This step fits a whole-genome regression model to account for population structure
# and relatedness by using all genotyped variants to predict the phenotype under
# the null hypothesis of no genetic associations.
#
# The output from this step (predicted values and model parameters) will be used
# in step 02 for association testing while accounting for this structure.
# --------------------------------

echo "[INFO] Start: $(date)"

check_dx_file() {
    # Checks the existence of a file
    # within DNAnexus.
    FILEPATH="$1"
    if ! dx ls "$FILEPATH" &> /dev/null; then
        echo "[ERROR] Missing file in DNAnexus: $FILEPATH $(date)"
        exit 1
    fi
}

download_dx_file() {
    FILEPATH="$1"
    OUTPATH="$2"
    check_dx_file "$FILEPATH"
    mkdir -p "$OUTPATH"
    dx download "$FILEPATH" -o "$OUTPATH" -f
}

# --------------------------------
# 00. Set up
# --------------------------------
set -euo pipefail
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <PHENOTYPE_CODE> <CHR (1-22)>"
    exit 1
fi

# Inputs
PHENO="$1"
CHR="$2"

if ! [[ "$CHR" =~ ^[0-9]+$ ]] || [[ "$CHR" -lt 1 ]] || [[ "$CHR" -gt 22 ]]; then
    echo "Error: Chromosome must be an integer between 1 and 22"
    exit 1
fi

CHR_PADDED=$(printf "%02d" "$CHR")

# --------------------------------
# 01. Check inputs and outputs
# MODIFY THESE PATHS TO MATCH YOUR PROJECT
# --------------------------------

INPUT_ROOT="/Users/Roberto/00_data"
GENO_DIR="${INPUT_ROOT}/01_genotypes/02_topmed-imputed/"
PHENO_DIR="${INPUT_ROOT}/00_phenotypes/02_mcps_qc-phenotypes/ver_2"
COVAR_DIR="${INPUT_ROOT}/00_phenotypes/02_mcps_qc-phenotypes/ver_2"
GRM_DIR="/Data/KING-IBD"

# Local
DOWNLOAD_DIR="${HOME}/downloads/regenie-inputs"
OUTPUT_ROOT="${HOME}/results"
OUTPUT_DIR="${OUTPUT_ROOT}/${PHENO}"

# File definitions
BGEN_FILE="${GENO_DIR}/op_prefix_chr${CHR}.shapeit5_ligated.high-quality.bgen"
PHENO_FILE="${PHENO_DIR}/${PHENO}.qc-phenotype.txt"
COVAR_FILE="${COVAR_DIR}/covars.qc-phenotype.nofasting.txt"
GRM_FILE="${GRM_DIR}/king_ibdseg_4th.seg"

# Downloading
echo "[INFO] $(date) | Downloading input files..."

download_dx_file "$PHENO_FILE" "$DOWNLOAD_DIR"
download_dx_file "$COVAR_FILE" "$DOWNLOAD_DIR"
download_dx_file "$GRM_FILE" "$DOWNLOAD_DIR"
download_dx_file "$BGEN_FILE" "$DOWNLOAD_DIR"

echo "[INFO] $(date) | All required inputs downloaded."

# Redefine the local inputs
BGEN_FILE="${DOWNLOAD_DIR}/op_prefix_chr${CHR}.shapeit5_ligated.high-quality.bgen"
PHENO_FILE="${DOWNLOAD_DIR}/${PHENO}.qc-phenotype.txt"
COVAR_FILE="${DOWNLOAD_DIR}/covars.qc-phenotype.nofasting.txt"
GRM_FILE="${DOWNLOAD_DIR}/king_ibdseg_4th.seg"

# Run the script
docker run -w /tmp/ \
	-v "${DOWNLOAD_DIR}":/input \
	-v "${OUTPUT_DIR}":/output \
	regenie:v3.0.1 \
	regenie \
		--step 1 \
		--bgen /input/$(basename "${BGEN_FILE}") \
		--covarFile /input/$(basename "${COVAR_FILE}") \
		--phenoFile /input/$(basename "${PHENO_FILE}") \
		--bt --bsize 1000 \
		--threads 4 \
		--gz --out /output/${PHENO}.gsav2_phased.unrelated \
		--verbose