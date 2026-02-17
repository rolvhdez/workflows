#!/bin/bash

# Author: Roberto Olvera Hernandez
# Date: 2026-02-16
# Description:
# The sample IDs of the phased VCF files of the MCPS genotype
# have 
#
# Usage:
# for i in {1..22}; do 
#   ./reformat_gsav2-chip.sh "/path/to/array/MCPS_Freeze_150.GT_hg38.pVCF.revised_qc_chr${i}_phased.vcf.gz"
# done

# --------------------------------
# 0. Set-up
# --------------------------------

set -euo pipefail
echo "[INFO] $(date) | Start"

# Inputs
vcf_file="$1" # input like *.phased.gz
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <input.vcf.gz>"
    exit 1
fi
prefix="${vcf_file%.vcf.gz}"
new_prefix="${prefix}.new-sample-ids"

# --------------------------------
# 1. Generate the mapping file
# --------------------------------

echo "[INFO] $(date) | Creating list of Sample IDs: "

sample_ids="/tmp/sample-ids.txt"
map_ids="/tmp/map-ids.txt"

# From first_file (arbitraty), get the sample IDs
bcftools query -l "$vcf_file" > "$sample_ids"
awk -F'-' 'BEGIN{OFS="\t"} {print $0, $2}' "$sample_ids" > "$map_ids"
head "$map_ids"
# if cut -f2 "$map_ids" | sort | uniq -d | grep -q .; then
#     echo "[ERROR] $(date) | Duplicate new sample IDs detected"
#     cut -f2 "$map_ids" | sort | grep "DUP"
#     exit 1
# fi
echo "[INFO] $(date) | Sample IDs: $(cat $map_ids | wc -l)"

# --------------------------------
# 2. Recode sample IDs
# --------------------------------

echo "[INFO] $(date) | Changing the headers"
new_vcf="/tmp/new-sample-ids.vcf.gz"
bcftools reheader -s "$map_ids" -o "$new_vcf" "$vcf_file"

# --------------------------------
# 3. Recode file and variant IDs
# --------------------------------

# NOTE:
# The flag
# --set-all-var-ids @:#:\$r:\$a
#

plink2 --vcf "$new_vcf" \
    --export bgen-1.2 \
    --snps-only \
    --rm-dup --geno 0.05 \
    --memory 3000 \
    --threads 1 \
    --set-all-var-ids @:#:\$r:\$a \
    --out "$new_prefix"

echo "[INFO] $(date) | Done!"