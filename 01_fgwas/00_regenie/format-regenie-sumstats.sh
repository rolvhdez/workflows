#!/bin/bash

# Author: Roberto Olvera Hernandez
# Date: 2026-03-06

# Description:
# This small script both combines the
# sumstats provided - outputs of this workflow -
# concatenates them into a single one,
# and then gives them a standard format.
#
# Usage:
# ls <your_files> | xargs ./format-regenie-sumstats.sh | gzip -c > <new_file>

set -euo pipefail

input_files=("$@")
first_file=true
temp_file=$(mktemp)
for file in "${input_files[@]}"; do
    if [[ $first_file == true ]]; then
        if [[ $file == *.gz ]]; then
            zcat "$file"
            first_file=false
        else
            cat "$file"
            first_file=false
        fi
    fi
    if [[ "$file" == *.gz ]]; then
        zcat "$file" | awk 'NR > 1'
    else
        awk 'NR > 1' "$file"
    fi
done |
    awk 'NR==1 {print $0, "P"; next} {print $0, 10^(-$13)}' | # Compute p-value
    awk 'NR==1 {print; next} {$3 = "chr"$3; print}' | # add 'chr' to SNP ID
    sed 's/ /\t/g' # change separator