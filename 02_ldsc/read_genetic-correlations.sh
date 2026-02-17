#!/bin/bash
set -e

out_dir="$HOME/results"
mkdir -p $out_dir
out_file="${out_dir}/rg_table.txt"

first_file=1
for file in $(ls $HOME/data/*.log | sort -V); do
    if [[ $first_file == 1 ]]; then
        cat $file | tail -n 5 | head -n 2 > $out_file
        first_file=0
    fi
        cat $file | tail -n 5 | head -n 2 | tail -n +2 >> $out_file
done