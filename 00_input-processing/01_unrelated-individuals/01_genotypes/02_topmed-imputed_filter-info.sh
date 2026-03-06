#!/bin/bash
# Author: Roberto Olvera Hernandez
# Date: 2026-03-04

set -euo pipefail

echo "[INFO] Start: $(date)"

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

# --------------------------------
# 02. Check input files
# --------------------------------

CHR="$1"
if ! [[ "$CHR" =~ ^[0-9]+$ ]] || [[ "$CHR" -lt 1 ]] || [[ "$CHR" -gt 22 ]]; then
    echo "Error: Chromosome must be an integer between 1 and 22"
    exit 1
fi
CHR_PADDED=$(printf "%02d" "$CHR")

BGEN="file-G36G1VQ0b67BZ5PK4ZFQ4qF9" # ID is provided below
SAMPLE="file-G46kY3805vBvbg8q6J5357B4" # ID is provided below
SNPLIST=$(dx find data --tag 'topmed-imputed' --tag 'snplist' --name '*.txt.gz') # INFO > 0.99

# Check inputs
inputs=("$BGEN" "$SAMPLE" "$SNPLIST")
echo "[INFO] $(date) | Checking input files..."
for i in "${inputs[@]}"; do dxCheckFile "$i"; done
echo "[INFO] $(date) | All required inputs exist."

# Create the output directory on DNAnexus
if ! dx ls "$OUT_DIR" &> /dev/null; then
    echo "[INFO] $(date) | Directory does not exist: $OUT_DIR"
    exit 1
fi

# --------------------------------
# 03. Launch Swiss Army Knife
# --------------------------------

JOB_NAME="Segment TOPMed-imputed (Unphased) INFO > 0.99 (chr_${CHR_PADDED})"
BGEN_NAME="$(dx describe $BGEN --json | jq -r .name)"
SAMPLE_NAME="$(dx describe $SAMPLE --json | jq -r .name)"
SNPLIST_NAME="$(dx describe $SNPLIST --json | jq -r .name)"

OUT_DIR="/Users/Roberto/00_data/01_genotypes/02_topmed-imputed/00_unphased/"
OUT_NAME="MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.high-quality.chr_${CHR}"

# DNAnexus Launch Command: Swiss Army Knife
# Instance: mem1_ssd1_v2_x8 (CPU : 8 | RAM : 16.0 GB | VOL : 200 GB)
dx run app-swiss-army-knife \
  -imount_inputs=true \
  -iin="$BGEN" \
  -iin="$SAMPLE" \
  -iin="$SNPLIST" \
  --instance-type mem1_ssd1_v2_x8 \
  --priority high \
  --name "$JOB_NAME" \
  --folder "$OUT_DIR" \
  -icmd="
  plink2 --bgen '${BGEN_NAME}' 'ref-first' \
    --sample '${SAMPLE_LIST}' \
    --extract '${SNPLIST_NAME}' \
    --chr '${CHR}' \
    --export bgen-1.2 \
    --memory 12000 --threads 8 \
    --out '${OUT_NAME}'
  " \
  --cost-limit 3 \
  --brief \
  --yes

# Retrieving the job name
echo "[DONE] $(date) | Job submitted: $JOB_NAME"
exit 0
