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

set -euo pipefail

# --------------------------------
# 0. Set-up
# --------------------------------
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <gsav2-chip.vcf.gz> <outdir/>"
    exit 1
fi

bed_prefix="$1"
output_dir="$2"

echo "[START] $(date) 'Reformating GSAv2-CHIP files'"

# --------------------------------
# 1. Convert to BED to BGEN with PLINK2
# --------------------------------
echo "[INFO] $(date) | Running PLINK2 conversion"
for i in {1..22}; do
    bed_name=$(basename $bed_prefix)
    output_name="${output_dir}/${bed_name}.chr${i}"
    plink2 --bfile "$bed_prefix" \
        --chr $i \
        --export bgen-1.2 \
        --rm-dup 'exclude-all' \
        --snps-only \
        --set-all-var-ids '@:#:$r:$a' \
        --out "${output_name}"
done
echo "[INFO] $(date) | Done!"