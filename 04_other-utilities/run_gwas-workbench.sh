#!/bin/bash
set -euo pipefail
shopt -s nullglob

cd $HOME/repositories/workflows/03_data-visualization/gwas-workbench

dict="$HOME/downloads/phenotype_dictionary.txt"

for file in "$HOME"/downloads/sumstats/*.gz; do
    filename=$(basename "$file")
    trait=${filename%%.*}
    description=$(grep $trait $dict | cut -d$'\t' -f2)
    #
    out_dir="$HOME/results/${trait}"
    mkdir -p $out_dir
    #
    Rscript gwas-make-plots.R "$file" "${out_dir}/${trait}" snipar no "$description" yes
done