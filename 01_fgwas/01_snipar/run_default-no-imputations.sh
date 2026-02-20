#!/bin/bash

# Author: Roberto Olvera Hernandez
# Date: 2026-02-20
#
# --------------------------------
# DEFAULT NO IMPUTATION
# --------------------------------
# Young et al. (2022) developed a model with 'snipar'
# to run a meta-analysis of a Within-Family
# GWAS model using sibling-pair and parent-offspring
# trios to obtain Direct Genetic Effects.
#
# To enhance power, they developed 'Mendelian Imputation',
# a method to obtain parental genotypes of missing
# parents in a cohort.
#
# This script launches a job to DNAnexus (Swiss Army Knife)
# for each chromosome for a given phenotype using
# 'snipar' with the sibpair/trio meta-analysis
# without providing imputed parental genotypes.
#
# https://doi.org/10.1038/s41588-025-02118-0
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

OUTPUT_ROOT="/Users/Roberto/01_results/01_gwas-sumstats/02_fgwas-snipar-topmed_phased/00_without-mendelian-imputation"
OUTPUT_DIR="${OUTPUT_ROOT}/batches/${PHENO}"

# File definitions
PHENO_FILE="${PHENO_DIR}/${PHENO}.qc-phenotype.txt"
COVAR_FILE="${COVAR_DIR}/covars.qc-phenotype.nofasting.txt"
GRM_FILE="${GRM_DIR}/king_ibdseg_4th.seg"

BGEN_FILE="${GENO_DIR}/op_prefix_chr${CHR}.shapeit5_ligated.high-quality.bgen"

# Check existence in DNAnexus
echo "[INFO] $(date) | Checking input files..."

check_dx_file "$PHENO_FILE"
check_dx_file "$COVAR_FILE"
check_dx_file "$GRM_FILE"
check_dx_file "$BGEN_FILE"

echo "[INFO] $(date) | All required inputs exist."

# Create the output directory on DNAnexus
if ! dx ls "$OUTPUT_DIR" &> /dev/null; then
    echo "[INFO] $(date) | Creating output directory: $OUTPUT_DIR"
    dx mkdir -p "$OUTPUT_DIR"
fi

# --------------------------------
# 02. Launch Swiss Army Knife
# --------------------------------

# Job name
JOB_NAME="FGWAS (default, no-imputation) ${PHENO} (chr_${CHR_PADDED})"

# DNAnexus Launch Command: Swiss Army Knife
dx run app-swiss-army-knife \
    -iimage="rolvhdez/snipar:0.0.22" \
    -imount_inputs=true \
    -iin="${PHENO_FILE}" \
    -iin="${COVAR_FILE}" \
    -iin="${GRM_FILE}" \
    -iin="${BGEN_FILE}" \
    --instance-type mem3_ssd1_v2_x32 \
    --priority high \
    --name "$JOB_NAME" \
    --folder "$OUTPUT_DIR" \
    -icmd="
    gwas.py '${PHENO}.qc-phenotype.txt' \
      --bgen 'op_prefix_chr@.shapeit5_ligated.high-quality' \
      --chr_range ${CHR} \
      --covar 'covars.qc-phenotype.nofasting.txt' \
      --grm 'king_ibdseg_4th.seg' \
      --sparse_thresh 0.075 \
      --min_maf 0.01 \
      --cpu 4 \
      --threads 1 \
      --batch_size 1000 \
      --out '${PHENO}.chr_@.topmed_phased.nonimputed'
    " \
    --brief
#    --yes

# Retrieving the job name
echo "[DONE] $(date) | Job submitted: $JOB_NAME"
exit 0
