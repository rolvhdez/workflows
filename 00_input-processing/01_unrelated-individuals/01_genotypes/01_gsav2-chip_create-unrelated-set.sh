#!/bin/bash
# Author: Roberto Olvera Hernandez (updated)
# Date: 2026-03-05
#
# Description:
# Reformats the BED with matching variant
# IDs with the TOPMed-Imputed (phased) dataset
# creating a new set of BGEN files.
# 
# Usage:
# ./reformat_gsav2-chip.sh <prefix_bfiles> <output_dir/>


set -euo pipefail

dxCheckFile() {
    # Checks the existence of a file
    # within DNAnexus.
    local FILE_ID="$1"
    if ! dx ls "$FILE_ID" &> /dev/null; then
        echo "[ERROR] Missing file in DNAnexus: $FILE_ID $(date)"
        exit 1
    fi
    echo "[INFO] Using: $FILE_ID"
}

echo "[INFO] Start: $(date)"

# --------------------------------
# 02. Check input files
# --------------------------------
BED="file-GzjfyXQ0fk00xkbVyy6J5f97"
BIM="file-GzjY8Zj0fk0B2gGkq8YPPqXg"
FAM="file-GzjY8900fk0PXYGp626pfj83"
UNRELATED=$(
dx find data \
  --tag 'unrelated-individuals' --tag 'sample-list' \
  --name "*.keep" --brief
)

# Check ---
echo "[INFO] $(date) | Checking input files..."
inputs=("$BED" "$BIM" "$FAM" "$UNRELATED")
for i in "${inputs[@]}"; do dxCheckFile "$i"; done
echo "[INFO] $(date) | All required inputs exist."

# --------------------------------
# 03. Launch Swiss Army Knife
# --------------------------------

JOB_NAME="Apply QC GSAv2-CHIP (Unrelated individuals) - Autosomes"

BED_NAME=$(dx describe "$BED" --json | jq -r .name)
BIM_NAME=$(dx describe "$BIM" --json | jq -r .name)
FAM_NAME=$(dx describe "$FAM" --json | jq -r .name)
UNRELATED_NAME=$(dx describe "$UNRELATED" --json | jq -r .name)

OUT_DIR="/Users/Roberto/00_data/01_genotypes/01_gsav2-chip/"
PREFIX="${BED_NAME%.bed}.maf-01.unrelated"

# Execute the launcher
dx run app-swiss-army-knife \
  -imount_inputs=true \
  -iin="$BED" \
  -iin="$BIM" \
  -iin="$FAM" \
  -iin="$UNRELATED" \
  --instance-type mem1_ssd1_v2_x4 \
  --priority high \
  --name "$JOB_NAME" \
  --folder "$OUT_DIR" \
  -icmd="
    for i in {1..22}; do
        outname="${PREFIX}\.chr\_\$i"
        plink2 --bfile "${BED_NAME%.bed}" \
            --chr \$i \
            --snps-only --rm-dup 'exclude-all' \
            --geno 0.05 --maf 0.01 \
            --set-all-var-ids 'chr@:#:\$r:\$a' \
            --export bgen-1.2 \
            --keep "${UNRELATED_NAME}" \
            --out "\$outname"
    done
  " \
  --cost-limit 3 \
  --brief \
  --yes