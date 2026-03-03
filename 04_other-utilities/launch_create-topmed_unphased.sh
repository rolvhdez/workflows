#!/bin/bash
# Author: Roberto Olvera Hernandez
# Date: 2026-03-03

set -euo pipefail

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
CHR="$1"
if ! [[ "$CHR" =~ ^[0-9]+$ ]] || [[ "$CHR" -lt 1 ]] || [[ "$CHR" -gt 22 ]]; then
    echo "Error: Chromosome must be an integer between 1 and 22"
    exit 1
fi
CHR_PADDED=$(printf "%02d" "$CHR")

# --------------------------------
# 02. Check input files
# --------------------------------

# Fixed files
SNP_LIST="file-J5k3G4j04QxY2gfZpqgPGzK4" # INFO > 0.99
SAMPLE_LIST="file-J5k3G4j04QxY2gfZpqgPGzK4" # Unrelated individuals

# Job name
BGEN="MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.bgen" # ID is provided below
SAMPLE="MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.COLLAB.sample" # ID is provided below
OUT_NAME="MCPS_Freeze_150.GT_hg38.pVCF.rgcpid.QC2.TOPMED_dosages.chr_${CHR}"
OUTPUT_ROOT="/Users/Roberto/00_data/01_genotypes/02_topmed-imputed/00_unphased/"

# --------------------------------
# 03. Launch Swiss Army Knife
# --------------------------------

JOB_NAME="Unrelated Individuals - TOPMed-imputed Unphased (chr_${CHR_PADDED})"

# DNAnexus Launch Command: Swiss Army Knife
dx run app-swiss-army-knife \
  -imount_inputs=true \
  -iin="${SNP_LIST}" \
  -iin="${SAMPLE_LIST}" \
  -iin="file-G36G1VQ0b67BZ5PK4ZFQ4qF9" \
  -iin="file-G46kY3805vBvbg8q6J5357B4" \
  --instance-type mem3_ssd2_x8 \
  --priority high \
  --name "$JOB_NAME" \
  --folder "$OUT_DIR" \
  -icmd="
  plink2 --bgen '${BGEN}' \
    --sample '${SAMPLE}' \
    --chr $CHR \
    --export bgen-1.2 \
    --keep '${SAMPLE_LIST}' \
    --extract '${SNP_LIST}' \
    --rm-dup 'exclude-all' \
    --snps-only \
    --set-all-var-ids 'chr\@:\#:\$r:\$a' \
    --out '${OUT_NAME}'
  "
