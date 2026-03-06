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
    FILE="$1"
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
PHENO_CODE="$1"
CHR="$2"

if ! [[ "$CHR" =~ ^[0-9]+$ ]] || [[ "$CHR" -lt 1 ]] || [[ "$CHR" -gt 22 ]]; then
    echo "Error: Chromosome must be an integer between 1 and 22"
    exit 1
fi

CHR_PADDED=$(printf "%02d" "$CHR")

# --------------------------------
# 01. Check input files
# --------------------------------

# File definitions
BGEN=$(dx find data --tag 'gsav2-chip' --tag 'unrelated-individuals' --name "*chr_${CHR}.bgen" --brief)
SAMPLE=$(dx find data --tag 'gsav2-chip' --tag 'unrelated-individuals' --name "*chr_${CHR}.sample" --brief)
PHENO=$(dx find data --tag 'phenotype' --tag 'unrelated-individuals' --name "${PHENO_CODE}*.txt" --brief)
COVAR=$(dx find data --tag 'covariates' --tag 'unrelated-individuals' --name "*.nofasting.txt" --brief)

# Check ---
inputs=("$BGEN" "$SAMPLE" "$PHENO" "$COVAR")
echo "[INFO] $(date) | Checking input files..."
for i in "${inputs[@]}"; do checkDxFile "$i"; done
echo "[INFO] $(date) | All required inputs exist."

# --------------------------------
# 03. Launch Swiss Army Knife
# --------------------------------

JOB_NAME="REGENIE (unrelated-individuals) ${PHENO_CODE} (chr_${CHR_PADDED})"

BGEN_NAME=$(dx describe $BGEN --json | jq -r .name)
SAMPLE_NAME=$(dx describe $SAMPLE --json | jq -r .name)
PHENO_NAME=$(dx describe $PHENO --json | jq -r .name)
COVAR_NAME=$(dx describe $COVAR --json | jq -r .name)

OUTPUT_ROOT="/Users/Roberto/01_results/01_gwas-sumstats/03_regenie-unrelated/batches"
OUTPUT_DIR="${OUTPUT_ROOT}/${PHENO_CODE}"

# Create the output directory on DNAnexus
if ! dx ls "$OUTPUT_DIR" &> /dev/null; then
    echo "[INFO] $(date) | Creating output directory: $OUTPUT_DIR"
    dx mkdir -p "$OUTPUT_DIR"
fi

# Run the script

dx run app-swiss-army-knife \
    -imount_inputs=true \
    -iin="${BGEN}" \
    -iin="${SAMPLE}" \
    -iin="${PHENO}" \
    -iin="${COVAR}" \
    --instance-type mem1_ssd2_v2_x16 \
    --priority high \
    --name "${JOB_NAME}" \
    --folder "${OUTPUT_DIR}" \
    -icmd="
	regenie \
		--step 1 \
		--bgen "${BGEN_NAME}" \
		--covarFile "${COVAR_NAME}" \
		--phenoFile "${PHENO_NAME}" \
		--bsize 1000 \
		--threads 4 \
		--gz --out "${PHENO_CODE}.gsav2-chip.unrelated.chr_${CHR}" \
		--verbose
    " \
    --cost-limit 3