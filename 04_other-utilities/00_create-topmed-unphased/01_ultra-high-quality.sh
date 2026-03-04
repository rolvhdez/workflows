#!/bin/bash
# Author: Roberto Olvera Hernandez
# Date: 2026-03-03

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

# Fixed files
BGEN="/Data/TOPMED-imputed/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.bgen" # ID is provided below
SAMPLE="/Data/TOPMED-imputed/MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.COLLAB.sample" # ID is provided below
SNPLIST="/Users/Roberto/00_data/01_genotypes/02_topmed-imputed/rgcpid.QC2.TOPMED_dosages.high-quality.snpstats.txt.gz" # INFO > 0.99

inputs=("$BGEN" "$SAMPLE" "$SNPLIST")

# Job name
OUT_DIR="/Users/Roberto/00_data/01_genotypes/02_topmed-imputed/00_unphased/"
OUT_NAME="MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.high-quality"

# Check inputs
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

JOB_NAME="Segment TOPMed-imputed (Unphased) INFO > 0.99"

# DNAnexus Launch Command: Swiss Army Knife
# Instance: mem1_ssd1_v2_x8 (CPU : 8 | RAM : 16.0 GB | VOL : 200 GB)
# Cost: _____ per hr
# Max runtime: _h (~$___)
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
  plink2 --bgen '$(basename $BGEN)' 'ref-first' \
    --extract '$(basename $SNPLIST)' \
    --export bgen-1.2 \
    --memory 12000 --threads 8 \
    --out '${OUT_NAME}'
  " \
  --cost-limit 3 \
  --brief

# Retrieving the job name
echo "[DONE] $(date) | Job submitted: $JOB_NAME"
exit 0
