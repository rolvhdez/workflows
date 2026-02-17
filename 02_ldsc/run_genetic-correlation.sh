#!/bin/bash

# Usage: ./script.sh <snipar_dir> <regenie_dir>

snipar_dir="$1"
regenie_dir="$2"

# Create unique temporary file names for each iteration to avoid conflicts
snipar_munged_base="/tmp/snipar_sumstats"
regenie_munged_base="/tmp/regenie_sumstats"

mkdir -p "$HOME/results/ldsc/"

# Function to extract phenotype ID from filename
get_phenotype_id() {
    local filename=$(basename "$1")
    # Remove everything after the first dot
    echo "${filename%%.*}"
}

munge_snipar() {
    local input_file="$1"
    local output_base="$2"
    local tmp_file="/tmp/snipar_${pheno_id}.temp.txt"
    
    echo "[DEBUG] Munging snipar file: $input_file"
    echo "[DEBUG] Output base: $output_base"
    echo "[DEBUG] Temp file: $tmp_file"
    
    # Create the P-value column
    zcat "$input_file" | awk 'BEGIN {OFS="\t"}
    {
        if (NR == 1) {
            # Add direct_P header
            print $0 "\tdirect_P"
        } else {
            # log10P should be column 11 (adjust if needed)
            log10P = $11
            # Convert -log10(P) to P
            P = 10^(-log10P)
            print $0 "\t" P
        }
    }' > "$tmp_file"
    
    echo "[DEBUG] Running munge_sumstats.py on temp file"
    # Check if temp file exists and has content
    if [ ! -s "$tmp_file" ]; then
        echo "[ERROR] Temp file $tmp_file is empty or doesn't exist"
        return 1
    fi
    
    # Run munge_sumstats
    $HOME/.local/bin/munge_sumstats.py \
        --sumstats "$tmp_file" \
        --snp "SNP" \
        --N-col "direct_N" \
        --signed-sumstats "direct_Beta,0" \
        --a1 "A1" \
        --a2 "A2" \
        --p "direct_P" \
        --out "$output_base"
    
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "[ERROR] munge_sumstats.py failed for snipar with exit code $exit_code"
    fi
    
    # Clean up temp file
    rm -f "$tmp_file"
    
    return $exit_code
}

# REGENIE ----
munge_regenie() {
    local input_file="$1"
    local output_base="$2"
    local tmp_file="${output_base}.tmp.sumstats"
    
    echo "[DEBUG] Munging regenie file: $input_file"
    echo "[DEBUG] Output base: $output_base"
    
    # Run munge_sumstats
    $HOME/.local/bin/munge_sumstats.py \
        --sumstats "$input_file" \
        --snp "ID" \
        --N-col "N" \
        --signed-sumstats "BETA,0" \
        --a1 "ALLELE1" \
        --a2 "ALLELE0" \
        --info "INFO" \
        --p "P" \
        --out "$output_base"
    
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "[ERROR] munge_sumstats.py failed for regenie with exit code $exit_code"
        return $exit_code
    fi
    
    # Check if output file was created
    if [ ! -f "${output_base}.sumstats.gz" ]; then
        echo "[ERROR] Output file ${output_base}.sumstats.gz not created"
        return 1
    fi
    
    # Add "chr" prefix to chromosome numbers
    echo "[DEBUG] Adding 'chr' prefix to SNP IDs"
    zcat "${output_base}.sumstats.gz" | \
        awk 'BEGIN {OFS="\t"} 
             NR==1 {print $0} 
             NR>1 && $1=="SNP" {print $0} 
             NR>1 && $1!="SNP" {print "chr"$0}' | \
        gzip > "${tmp_file}.gz"
    
    # Move the modified file back
    mv "${tmp_file}.gz" "${output_base}.sumstats.gz"
    
    return 0
}

# Heritability analysis ---
run_ldsc_rg() {
    echo "[DEBUG] Running ldsc.py with files:"
    echo "[DEBUG]   File 1: $1"
    echo "[DEBUG]   File 2: $2"
    echo "[DEBUG]   LD scores: $3"
    echo "[DEBUG]   Output: $4"
    
    # Check if input files exist
    if [ ! -f "$1" ]; then
        echo "[ERROR] Input file 1 not found: $1"
        return 1
    fi
    if [ ! -f "$2" ]; then
        echo "[ERROR] Input file 2 not found: $2"
        return 1
    fi
    
    $HOME/.local/bin/ldsc.py \
        --rg "$1,$2" \
        --ref-ld-chr "$3" \
        --w-ld-chr "$3" \
        --out "$4"
    
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "[ERROR] ldsc.py failed with exit code $exit_code"
    fi
    
    return $exit_code
}

# Process all snipar sumstats files
for snipar_file in "$snipar_dir"/*.sumstats.gz; do
    # Skip if no files found
    if [ ! -f "$snipar_file" ]; then
        echo "[WARNING] No snipar files found in $snipar_dir"
        continue
    fi
    
    # Get the base filename without extension
    pheno_id=$(get_phenotype_id "$snipar_file")
    echo "[INFO] Processing phenotype: $pheno_id"
    
    # Create unique temporary file names for this iteration
    snipar_munged="${snipar_munged_base}_${pheno_id}"
    regenie_munged="${regenie_munged_base}_${pheno_id}"
    
    # Clean up any existing temporary files
    rm -f "${snipar_munged}.sumstats.gz" "${regenie_munged}.sumstats.gz"
    
    # Find corresponding regenie file
    # This looks for any regenie file that starts with the same phenotype ID
    regenie_file=$(find "$regenie_dir" -name "${pheno_id}.*.sumstats.gz" -type f 2>/dev/null | head -1)
    
    if [ -f "$regenie_file" ]; then
        echo "[INFO] Processing pair:"
        echo "  Snipar: $snipar_file"
        echo "  Regenie: $regenie_file"
        echo "  Phenotype ID: $pheno_id"
        
        # Set output name
        outname="$HOME/results/ldsc/${pheno_id}.genetic-correlations"
        
        # 1. Reformat the summary statistics
        echo "[INFO] Reformating $snipar_file"
        if ! munge_snipar "$snipar_file" "$snipar_munged"; then
            echo "[ERROR] Failed to munge snipar file: $snipar_file"
            continue
        fi
        
        echo "[INFO] Reformating $regenie_file"
        if ! munge_regenie "$regenie_file" "$regenie_munged"; then
            echo "[ERROR] Failed to munge regenie file: $regenie_file"
            # Clean up snipar file if regenie failed
            rm -f "${snipar_munged}.sumstats.gz"
            continue
        fi
        
        # Check if munged files were created
        if [ ! -f "${snipar_munged}.sumstats.gz" ]; then
            echo "[ERROR] Munged snipar file not created: ${snipar_munged}.sumstats.gz"
            rm -f "${regenie_munged}.sumstats.gz"
            continue
        fi
        
        if [ ! -f "${regenie_munged}.sumstats.gz" ]; then
            echo "[ERROR] Munged regenie file not created: ${regenie_munged}.sumstats.gz"
            rm -f "${snipar_munged}.sumstats.gz"
            continue
        fi
        
        # 2. Run the heritability analysis
        echo "[INFO] Writing to $outname"
        if ! run_ldsc_rg "${snipar_munged}.sumstats.gz" "${regenie_munged}.sumstats.gz" "$HOME/data/ldscores/MCPS_WGS_covLD_scores_1Mb_chr@" "$outname"; then
            echo "[ERROR] LDSC analysis failed for $pheno_id"
        fi
        
        # Clean up temporary files
        rm -f "${snipar_munged}.sumstats.gz" "${regenie_munged}.sumstats.gz"
        
        echo "[INFO] Completed processing for $pheno_id"
        echo "----------------------------------------"
    else
        echo "[WARNING] No matching regenie file found for $snipar_file"
        echo "[WARNING] Searched for: $regenie_dir/${pheno_id}.*.sumstats.gz"
        # List what files are actually in the regenie directory
        echo "[DEBUG] Files in $regenie_dir:"
        ls -1 "$regenie_dir"/*.sumstats.gz 2>/dev/null | head -5 || echo "  (no .sumstats.gz files found)"
    fi
done

echo "[INFO] Processing complete!"