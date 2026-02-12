#!/bin/bash

DIR="$1"

table_h2_logs() {
    local pattern=$1
    local log_files=($pattern) # expand the regex

    headers="trait\th2\th2_se\tlambda\tchi\tintercept\tintercept_se"
    #
    echo -e "$headers"
    for file in ${log_files[@]}; do
        # 1. name of trait
        trait=$(basename $file)
        trait="${trait%.log}"

        # 2. rest of fields
        h2=$(grep -E 'Total Observed scale h2:' "$file" | grep -oE '[0-9]+\.[0-9]+' | head -1)
        h2_se=$(grep -E 'Total Observed scale h2:' "$file" | grep -oE '[0-9]+\.[0-9]+' | tail -n 1)
        inflation=$(grep -E 'Lambda GC:' "$file" | grep -oE '[0-9]+\.[0-9]+')
        chi=$(grep -E 'Mean Chi\^2:' "$file" | grep -oE '[0-9]+\.[0-9]+')
        intercept=$(grep -E 'Intercept:' "$file" | grep -oE '[0-9]+\.[0-9]+' | head -1)
        intercept_se=$(grep -E 'Intercept:' "$file" | grep -oE '[0-9]+\.[0-9]+' | tail -1)

        # 3. make a row
        echo -e "$trait\t$h2\t$h2_se\t$inflation\t$chi\t$intercept\t$intercept_se"
    done
}

table_h2_logs "$DIR/*.log"