# Author: Roberto Olvera Hernandez
# Date: 2025-11-07
#
# Description:
#
# Note:
# Make sure that you're using your inputs IN ORDER.
#
# Usage:
# R does not support (at least not easy) flagged
# arguments for parsing. So, this is the equivalent
# of a `-h` of the script:
#
# Supported summary statistics formats:
# - snipar (v0.0.22) 
# - REGENIE
#
# --sumstats:path
#       Path to the summary statistics (.gz)
#
# --output:path
#       Output directory. If the directory does not exist,
#       the script will make a "prefix" using the name of the directory
#       you're trying to submit.
#       (Example) For `$HOME/path/001` the outputs will be `$HOME/path/001.qqplot.png`
#
# --model:str
#       Format of the summary statistics.
#       Supported: "snipar", "regenie"
#       Default: "snipar"
#
# --compute_bonferroni:str (optional)
#       If the script should compute the Bonferroni correction 
#       given the number of variants in the summary statistics.
#       Default: "no" (ie. p > 5e-8)
#
# --phenotype:str (optional)
#       Name of the phenotype to be used as title of the plots.
#       If not provided, the plots will not have a title.
#
# --annotations:path (optional)
#       Will perform annotations for the plots.
#       If not provided, the plots will not have annotations.
#       Default: "no".

suppressPackageStartupMessages(library(cli))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

# Load personalized functions
"%&%" <- function(a, b) paste0(a, b)
source("utils/00_utils.R")
source("utils/01_map_genes.R")
source("utils/02_qq_plot.R")
source("utils/03_manhattan_plot.R")
source("utils/04_effect_sizes.R")

# Process arguments
args <- commandArgs(trailingOnly = TRUE)

sumstats_file <- args[1]
output_dir <- args[2]
model <- if (length(args) < 3) "snipar" else args[3] # Default: snipar (v0.0.22)
compute_bonferroni <- if (length(args) < 4) "no" else args[4] # Default: empty
phenotype <- if (length(args) < 5) "" else args[5] # Default: empty
annotations <- if (length(args) < 6) "no" else args[6] # Default: empty

# Parameter error handling
if (!file.exists(sumstats_file)) {
  # Check if the summary statistics exist
  cli_abort(c(
    "{sumstats_file} does not exist",
    "x" = "You've supplied a file that does not exist."
  ))
}
if (!dir.exists(output_dir)) {
  # Provide a file in the directory for output
  output_dir <- output_dir %&% "."
}

#------------------------------------
# Read the summary statistics ---
raw_sumstats <- fancy_process(
  process = read_sumstats_file,
  message = "Reading " %&% sumstats_file,
  # Function parameters
  sumstats_path = sumstats_file,
  chunk_size = 1000000
)
df_sumstats <- reformat_sumstats(raw_sumstats, model) # Reformat the table
df_sumstats$CHR <- as.integer(df_sumstats$CHR)
k <- length(unique(df_sumstats$SNP)) # Number of lines
if (compute_bonferroni == "yes") {
  bonferroni <- 0.05 / nrow(df_sumstats)
} else {
  bonferroni <- 5e-8 # Bonferroni adjusted P-Value
}
cli_alert_info(scales::comma(k) %&% " SNPs found in `" %&% sumstats_file %&% "`.")
cli_alert_info("Bonferroni adjusted P-value: " %&% scales::scientific(bonferroni))

# Make the annotations ---
if (annotations == "yes") {
  # Check that there are significant SNPs to annotate
  sig_k <- df_sumstats %>% filter(P <= bonferroni) %>% pull(SNP)
  if (length(sig_k) > 0) {
    cli::cli_alert_warning(scales::comma(length(sig_k)) %&% " SNPs found at p <=" %&% bonferroni)
    genes <- create_gene_ranges()
    df_annotations <- annotate_genes_to_sig_snps(
      sumstats = df_sumstats,
      gene_range = genes,
      sig = bonferroni
    )
    df_annotations <- df_annotations %>% select(SNP, GENE)
    df_sumstats <- df_sumstats %>%
      left_join(df_annotations, by = "SNP")
  } else {
    cli::cli_alert_warning("No significant SNPs were found at p <= " %&% bonferroni %&% ". Skipping annotation.")
  }
}
print(head(df_sumstats))

### QQ PLOT ### ------------------------------------------
qq_list <- get_qqvalues(df_sumstats)
pvalues <- qq_list[[1]]
lambda <- qq_list[[2]]
qq_plot <- make_qqplot(pvalues, phenotype, lambda)
cli_alert_info("Lambda genetic inflation factor: " %&% round(lambda, 4))
export_plot(qq_plot, output_dir %&% "qqplot.png")

### MANHATTAN PLOT ### ------------------------------------------
list_manhattan <- format_for_manhattan(df_sumstats)
df_manhattan <- list_manhattan[[1]]
df_axis <- list_manhattan[[2]]
manhattan_plot <- make_manhattan(df_manhattan, df_axis, phenotype, bonferroni)
if ("GENE" %in% names(df_sumstats)) {
  # Add the annotation labels
  manhattan_plot <- manhattan_plot +
    geom_text_repel(
      data = subset(df_manhattan, !is.na(df_manhattan$GENE)),
      aes(label = GENE),
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = 16,
      size = 3
    )
}
export_plot(manhattan_plot, output_dir %&% "manhattan_plot.png")

### EFFECT SIZES ### ------------------------------------------
# Effect sizes vs. p-values
effect_pvalue_plot <- make_effectsizes_plot(df_sumstats, phenotype, bonferroni)
if ("GENE" %in% names(df_sumstats)) {
  # Add the annotation labels
  effect_pvalue_plot <- effect_pvalue_plot +
    geom_text_repel(
      data = subset(df_sumstats, !is.na(df_sumstats$GENE)),
      aes(
        x = BETA,
        y = -log10(P),
        label = GENE
      ),
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = 16,
      size = 3
    )
}
export_plot(effect_pvalue_plot, output_dir %&% "effects_pvalues.png")

# Effect sizes vs. MAF
if ("MAF" %in% names(df_sumstats)) {
  effects_maf_plot <- make_effectmaf_plot(df_sumstats, phenotype, bonferroni)
  if ("GENE" %in% names(df_sumstats)) {
    # Add the annotation labels
    effects_maf_plot <- effects_maf_plot +
      geom_text_repel(
        data = subset(df_sumstats, !is.na(df_sumstats$GENE)),
        aes(
          x = MAF,
          y = BETA,
          label = GENE
        ),
        box.padding = 0.5,
        point.padding = 0.3,
        max.overlaps = 16,
        size = 3
      )
  }
  export_plot(effects_maf_plot, output_dir %&% "effects_maf.png")
}