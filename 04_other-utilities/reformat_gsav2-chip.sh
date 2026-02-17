#!/bin/bash
# Author: Roberto Olvera Hernandez (updated)
# Date: 2026-02-16
# Description:
#   Renames sample IDs in a phased VCF from "A-A" to "A".
# Usage:
#   for i in {1..22}; do
#     ./reformat_gsav2-chip.sh "/path/to/array/MCPS_Freeze_150.GT_hg38.pVCF.revised_qc_chr${i}_phased.vcf.gz"
#   done

set -euo pipefail

# --------------------------------
# 0. Set-up
# --------------------------------
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <input.vcf.gz>"
    exit 1
fi

vcf_file="$1"
prefix="${vcf_file%.vcf.gz}"
new_prefix="${prefix}.new-sample-ids"

echo "[INFO] $(date) | Start processing $vcf_file"

# Create a temporary directory that will be cleaned up on exit
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

# --------------------------------
# 1. Generate and validate mapping
# --------------------------------
echo "[INFO] $(date) | Extracting sample IDs"

sample_ids="$tmpdir/sample-ids.txt"
map_ids="$tmpdir/map-ids.txt"

bcftools query -l "$vcf_file" > "$sample_ids"
n_total=$(wc -l < "$sample_ids")
echo "[INFO] $(date) | Total samples in VCF: $n_total"

# Check that every sample ID contains a hyphen
if grep -qv '-' "$sample_ids"; then
    echo "[ERROR] $(date) | Some sample IDs lack a hyphen:"
    grep -v '-' "$sample_ids" | head -5
    exit 1
fi

# Build mapping: original ID -> second field after '-'
awk -F'-' 'BEGIN{OFS="\t"} {print $0, $2}' "$sample_ids" > "$map_ids"

# Verify that the new IDs are all non‑empty (should be guaranteed by hyphen check)
if ! awk '$2 == "" {c++} END{exit c}' "$map_ids"; then
    echo "[ERROR] $(date) | Some new IDs are empty (should not happen)"
    exit 1
fi

# Check for duplicate new IDs
duplicates=$(cut -f2 "$map_ids" | sort | uniq -d)
if [[ -n "$duplicates" ]]; then
    echo "[ERROR] $(date) | Duplicate new sample IDs detected:"
    echo "$duplicates"
    exit 1
fi

echo "[INFO] $(date) | Mapping valid for $n_total samples"

# --------------------------------
# 2. Recode sample IDs with bcftools
# --------------------------------
echo "[INFO] $(date) | Renaming samples with bcftools reheader"
new_vcf="$tmpdir/new-sample-ids.vcf.gz"
bcftools reheader -s "$map_ids" -o "$new_vcf" "$vcf_file"

# --------------------------------
# 3. Convert to BGEN with PLINK2
# --------------------------------
echo "[INFO] $(date) | Running PLINK2 conversion"
plink2 --vcf "$new_vcf" \
    --export bgen-1.2 \
    --snps-only \
    --rm-dup \
    --geno 0.05 \
    --memory 3000 \
    --threads 1 \
    --set-all-var-ids '@:#:\$r:\$a' \
    --out "$new_prefix"

echo "[INFO] $(date) | Done!"