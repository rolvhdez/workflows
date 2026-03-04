#!/bin/bash
# Author: Roberto Olvera Hernandez
# Date: 2026-03-04

set -euo pipefail

echo "[INFO] Start: $(date)"

dxCheckFile() {
    # Checks the existence of a file
    # within DNAnexus.
    FILEPATH="$1"
    if ! dx ls "$FILEPATH" &> /dev/null; then
        echo "[ERROR] Missing file in DNAnexus: $FILEPATH $(date)"
        exit 1
    fi
    echo "[INFO] $(date) | Ready: $(basename $FILEPATH) "
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

GENO_DIR="/Users/Roberto/00_data/01_genotypes/02_topmed-imputed/00_unphased"
BGEN="${GENO_DIR}/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.high-quality.bgen"
SAMPLE="${GENO_DIR}/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.high-quality.sample"
inputs=("$BGEN" "$SAMPLE")

OUT_DIR="/Users/Roberto/00_data/01_genotypes/02_topmed-imputed/00_unphased/"
OUT_NAME="MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.high-quality.applied-qc.chr_${CHR}"

# Check ---
echo "[INFO] $(date) | Checking input files..."
for i in "${inputs[@]}"; do dxCheckFile "$i"; done
echo "[INFO] $(date) | All required inputs exist."

if ! dx ls "$OUT_DIR" &> /dev/null; then
    echo "[INFO] $(date) | Directory does not exist: $OUT_DIR"
    exit 1
fi

# --------------------------------
# 03. Launch Swiss Army Knife
# --------------------------------

JOB_NAME="Apply QC TOPMed-imputed (Unphased) (chr_${CHR_PADDED})"

dx run app-swiss-army-knife \
  -imount_inputs=true \
  -iin="$BGEN" \
  -iin="$SAMPLE" \
  -iin="$CHR" \
  --instance-type mem1_ssd1_v2_x8 \
  --priority high \
  --name "$JOB_NAME" \
  --folder "$OUT_DIR" \
  -icmd="
  plink2 --bgen '$(basename $BGEN)' 'ref-first' \
    --sample '$(basename $SAMPLE)' \
    --chr '${CHR}'
    --snps-only --rm-dup 'exclude-all' \
    --geno 0.05 --maf 0.01 \
    --export bgen-1.2 \
    --memory 12000 --threads 8 \
    --out '${OUT_NAME}'
  " \
  --cost-limit 3 \
  --brief
