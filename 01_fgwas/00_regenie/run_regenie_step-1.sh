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

checkDxFile() {
    # Checks the existence of a file
    # within DNAnexus.
    FILEPATH="$1"
    if ! dx ls "$FILE" &> /dev/null; then
        echo "[ERROR] $(date) | Missing file in DNAnexus: $FILEPATH"
        exit 1
    fi
    echo "[INFO] Using: $FILE"
}
downloadDxFile() {
    FILE="$1"
    OUTPATH="$2"
    checkDxFile "$FILE"
    mkdir -p "$OUTPATH"
    dx download "$FILE" -o "$OUTPATH" -f
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

# Local
DOWNLOAD_DIR="${HOME}/downloads/regenie-inputs"
OUTPUT_ROOT="${HOME}/results"
OUTPUT_DIR="${OUTPUT_ROOT}/${PHENO}"

# File definitions
BGEN=$(dx find data --tag 'gsav2-chip' --tag 'unrelated-individuals' --name "*chr_${CHR}.bgen" --brief)
SAMPLE=$(dx find data --tag 'gsav2-chip' --tag 'unrelated-individuals' --name "*chr_${CHR}.sample" --brief)
PHENO=$(dx find data --tag 'phenotype' --tag 'unrelated-individuals' --name "${PHENO}*.txt" --brief)
COVAR=$(dx find data --tag 'covariates' --tag 'unrelated-individuals' --name "*.nofasting.txt" --brief)

# Downloading
echo "[INFO] $(date) | Downloading input files..."

downloadDxFile "$BGEN" "$DOWNLOAD_DIR"
downloadDxFile "$SAMPLE" "$DOWNLOAD_DIR"
downloadDxFile "$PHENO" "$DOWNLOAD_DIR"
downloadDxFile "$COVAR" "$DOWNLOAD_DIR"

echo "[INFO] $(date) | All required inputs downloaded."

# Redefine the local inputs
BGEN_NAME=$(dx describe $BGEN --json | jq -r .name)
SAMPLE_NAME=$(dx describe $SAMPLE --json | jq -r .name)
PHENO_NAME=$(dx describe $PHENO --json | jq -r .name)
COVAR_NAME=$(dx describe $COVAR --json | jq -r .name)

# Run the script
docker run -w /tmp/ \
	-v "${DOWNLOAD_DIR}":/input \
	-v "${OUTPUT_DIR}":/output \
	ghcr.io/rgcgithub/regenie/regenie:v3.0.1.gz \
	regenie \
		--step 1 \
		--bgen "/input/${BGEN_NAME}" \
		--covarFile "/input/${COVAR_NAME}" \
		--phenoFile "/input/${PHENO_NAME}" \
		--bsize 1000 \
		--threads 4 \
		--gz --out /output/${PHENO}.chr_${CHR}.gsav2_phased.unrelated \
		--verbose