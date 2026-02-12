#!/bin/bash
set -e

# Author: Roberto Olvera Hernandez
# Date: 2025-09-24
#
# Description:
# This small script allows to merge row-by-row a 
# selection of text files (compressed or not) into a single 
# 'uncompressed' text file.
#
# Usage (E01):
# ./merge-files.sh file_chr1.txt.gz file_chr2.txt.gz ... file_chr22.txt.gz > new_file.txt && gzip new_file.txt
#
# Usage (E02): 
# ./merge-files.sh file_chr*.txt.gz > new_file.txt && gzip new_file.txt
#
# Note:
# The script will assume the first file it reads contains a header
# and it will use the same one for the rest of the files.

input_files=("$@")
first_file=true
for file in ${input_files[@]}; do
    if [[ $first_file == true ]]; then
        if [[ $file == *.gz ]]; then
            zcat "$file" | head -1
            first_file=false
        else
            head "$file"
            first_file=false
        fi
    fi
    if [[ "$file" == *.gz ]]; then
        zcat "$file" | awk 'NR > 1'
    else
        awk 'NR > 1' "$file"
    fi
done