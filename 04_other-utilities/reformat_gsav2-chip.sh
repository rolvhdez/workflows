#!/bin/bash
# Author: Roberto Olvera Hernandez (updated)
# Date: 2026-02-17
#
# Description:
# Reformats the BED with matching variant
# IDs with the TOPMed-Imputed (phased) dataset
# creating a new set of BGEN files.
# 
# Usage:
# ./reformat_gsav2-chip.sh <prefix_bfiles> <output_dir/>

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
# 0. Set-up
# --------------------------------
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <gsav2-chip.vcf.gz>"
    exit 1
fi

# Inputs
bed_prefix="$1"

# --------------------------------
# 01. Check inputs and outputs
# MODIFY THESE PATHS TO MATCH YOUR PROJECT
# --------------------------------

INPUT_ROOT="/Data"
GENO_DIR="${INPUT_ROOT}/GSAv2-Chip/data/pVCF/revised_qc"
BED_FILE="${GENO_DIR}/${bed_prefix}.bed"
BIM_FILE="${GENO_DIR}/${bed_prefix}.bim"
FAM_FILE="${GENO_DIR}/${bed_prefix}.fam"

OUTPUT_ROOT="/Users/Roberto/00_data/01_genotypes"
OUTPUT_DIR="${OUTPUT_ROOT}/01_gsav2-chip"

# Check existence in DNAnexus
echo "[INFO] $(date) | Checking input files..."

check_dx_file "$BED_FILE"
check_dx_file "$BIM_FILE"
check_dx_file "$FAM_FILE"

# Create the output directory on DNAnexus
if ! dx ls "$OUTPUT_DIR" &> /dev/null; then
    echo "[INFO] $(date) | Creating output directory: $OUTPUT_DIR"
    dx mkdir -p "$OUTPUT_DIR"
fi

# --------------------------------
# 02. Launch Swiss Army Knife
# --------------------------------

# Job name
JOB_NAME="BED to BGEN - Homogenize variant ID's and minimal QC"

# DNAnexus Launch Command: Swiss Army Knife
dx run app-swiss-army-knife \
    -imount_inputs=true \
    -iin="${BED_FILE}" \
    -iin="${BIM_FILE}" \
    -iin="${FAM_FILE}" \
    --instance-type mem1_ssd2_v2_x8 \
    --priority high \
    --name "$JOB_NAME" \
    --folder "$OUTPUT_DIR" \
    -icmd="
    # Get the base name without extension
    base_name=\$(basename ${BED_FILE} .bed)
    
    for i in {1..22}; do
        plink2 --bfile \$base_name \
            --chr \$i \
            --export bgen-1.2 \
            --rm-dup 'exclude-all' \
            --snps-only \
            --geno 0.05 --maf 0.1 \
            --set-all-var-ids 'chr@:#:\$r:\$a' \
            --out \${base_name}-maf01.chr_\$i
    done
    " \
    --brief \
    --yes

# Retrieving the job name
echo "[DONE] $(date) | Job submitted: $JOB_NAME"