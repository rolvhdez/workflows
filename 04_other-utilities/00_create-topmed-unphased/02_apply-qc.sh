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

# Find the ID's of the required files
BGEN=$(
dx find data \
  --tag 'topmed-imputed' --tag 'unphased' --tag 'high-quality' \
  --name "*chr_${CHR}.bgen" --brief
) 
SAMPLE=$(
dx find data \
  --tag 'topmed-imputed' --tag 'unphased' --tag 'high-quality' \
  --name "*chr_${CHR}.sample" --brief
) 
UNRELATED=$(
dx find data \
  --tag 'unrelated-individuals' --tag 'sample-list' \
  --name "*.keep" --brief
)

# Check ---
echo "[INFO] $(date) | Checking input files..."

inputs=("$BGEN" "$SAMPLE" "$UNRELATED")
for i in "${inputs[@]}"; do dxCheckFile "$i"; done
echo "[INFO] $(date) | All required inputs exist."

# --------------------------------
# 03. Launch Swiss Army Knife
# --------------------------------

JOB_NAME="Apply QC TOPMed-imputed (Unphased) (chr_${CHR_PADDED})"

BGEN_NAME=$(dx describe "$BGEN" --json | jq -r .name)
SAMPLE_NAME=$(dx describe "$SAMPLE" --json | jq -r .name)
UNRELATED_NAME=$(dx describe "$UNRELATED" --json | jq -r .name)

# Output directory and file name
OUT_DIR="/Users/Roberto/00_data/01_genotypes/02_topmed-imputed/00_unphased/"
OUT_NAME="MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.high-quality.maf-01.unrelated.chr_${CHR}"
if ! dx ls "$OUT_DIR" &> /dev/null; then
    echo "[INFO] $(date) | Directory does not exist: $OUT_DIR"
    exit 1
fi

# Execute the launcher
dx run app-swiss-army-knife \
  -imount_inputs=true \
  -iin="$BGEN" \
  -iin="$SAMPLE" \
  -iin="$UNRELATED" \
  --instance-type mem1_ssd1_v2_x8 \
  --priority high \
  --name "$JOB_NAME" \
  --folder "$OUT_DIR" \
  -icmd="
  #
  plink2 --bgen '${BGEN_NAME}' 'ref-first' \
    --sample '${SAMPLE_NAME}' \
    --keep '${UNRELATED_NAME}' \
    --snps-only --rm-dup 'exclude-all' \
    --geno 0.05 --maf 0.01 \
    --make-pgen \
    --out 'tmpfile' && \
  plink2 --pfile 'tmpfile' \
    --export bgen-1.2 \
    --out '${OUT_NAME}'
  " \
  --cost-limit 3 \
  --yes \
  --brief
