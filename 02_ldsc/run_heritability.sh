#!/bin/bash
set -e

# Reformat summary statistics ----
# SNIPAR ----
munge_snipar() {
$HOME/.local/bin/munge_sumstats.py \
    --sumstats $1 \
    --snp "SNP" \
    --N-col "direct_N" \
    --signed-sumstats "direct_Beta,0" \
    --a1 "A1" \
    --a2 "A2" \
    --p "direct_P" \
    --out $2
}
# REGENIE ----
munge_regenie() {
    $HOME/.local/bin/munge_sumstats.py \
        --sumstats $1 \
        --snp "ID" \
        --N-col "N" \
        --signed-sumstats "BETA,0" \
        --a1 "ALLELE1" \
        --a2 "ALLELE0" \
        --info "INFO" \
        --p "P" \
        --out $2
}
# Heritability analysis ---
run_ldsc_h2() {
    $HOME/.local/bin/ldsc.py \
        --h2 $1 \
        --ref-ld-chr $2 \
        --w-ld-chr $2 \
        --out $3
}

# Main ----
sumstats_path="$1"
model="$2"
ldscores="$3"
munged_sumstats="/tmp/munged_sumstats"

mkdir -p "$HOME/results/ldsc/"
outname="$(basename $sumstats_path)"
outname="$HOME/results/ldsc/${outname%.gz}"

# 1. Reformat the summary statistics ---
echo "[INFO] Reformating $sumstats_path"

if [[ $model == "snipar" ]]; then
    # Add P value column
    tmp_file="/tmp/snipar.temp.txt"
    zcat $1 | awk 'BEGIN {OFS=" "}
    {
        if (NR == 1) {
            print $0 " direct_P"
        } else {
            log10P = $11
            P = 10^(-log10P)
            print $0 " " P
        }
    }' > $tmp_file
    #
    munge_snipar $tmp_file $munged_sumstats
    rm -f $tmp_file
elif [[ $model == "regenie" ]]; then
    munge_regenie $1 $munged_sumstats

    echo "[INFO] Adding 'chr' prefix to SNP column for 'regenie' format"
    tmp_file="${munged_sumstats}.tmp.sumstats"
    zcat "${munged_sumstats}.sumstats.gz" | \
        awk 'BEGIN {OFS="\t"} {if (NR==1 && $1=="SNP") print $0; else if (NR>1) print "chr"$0; else print $0}' | \
        gzip > "${tmp_file}.gz"
    mv "${tmp_file}.gz" "${munged_sumstats}.sumstats.gz"
    rm -f "$tmp_file"
else
    echo "[ERROR] Unknown model $model. Use 'snipar' or 'regenie'." >&2
    exit 1
fi

# 2. Run the heritability analysis
echo "[INFO] Writing to $outname"
run_ldsc_h2 "${munged_sumstats}.sumstats.gz" "$ldscores" "$outname"
rm "${munged_sumstats}.sumstats.gz"