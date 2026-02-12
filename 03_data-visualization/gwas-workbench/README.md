# Analysis Room

- **Author:** Roberto Olvera Hernandez
- **Date:** 2025-11-07
- **Based on:** https://github.com/rolvhdez/analysis-room

## Overview

This R scripts generate visualization for GWAS summary statistics and optionally performs variant annotation and Bonferroni correction.

It supports multiple summary statistics formats and is designed for command-line execution.

> Important:
> R does not natively support standard flagged argument parsing in a straightforward way.
> Therefore, all inputs must be provided in the exact order expected by the script.

Current supported formats:

- snipar (v0.0.22)
- REGENIE

## Usage
### 1. Manhattan Plots and more (`gwas-make-plots.R`)

Arguments are expected in order:

1. `sumstats` Path to the compressed summary statistics file (`.gz`).
2. `output` Output directory. If the directory does not exist, the script will generate a prefix based on the provided directory name - avoid using a `/` if you want the prefix.
3. `model` Format of the summary statistics file.
4. `compute_bonferroni` Wether to compute a Bonferroni correction or not; use genome-wide significance (p < 5e-8). Default is not.
5. `phenotype` Name of the phenotype to use. Default is `NULL`.
6. `annotations` Wether to annotate significant variants or not. It downloads ENSEMBL's GRCh38 version.

```shell
# Example (example: 001.manhattan.png)
Rscript gwas-make-plots.R /path/to/sumstats.gz /path/to/results/001 regenie no "Height" yes
```

### 2. Effect size comparisons (`effect-size-comparison.R`)

Arguments are presented in order:

1. `sumstats_a/b` Path to compressed summary statistics. The first one will be the one on the X axis of the scatter plot.
2. `output` Output directory. If the directory does not exist, the script will generate a prefix based on the provided directory name - avoid using a `/` if you want the prefix.
3. `phenotype` Name of the phenotype to use. Default is `NULL`.
4. `model_a/b` Format of the summary statistics files.
5. `model_label_a/b` Labels to be used to identify the models.

```shell
# Example (using snipar and regenie)
Rscript effect-size-comparison.R /path/to/regenie-sumstats.gz /path/to/snipar-sumstats.gz /path/to/results/001 "Height" regenie snipar "REGENIE" "SNIPAR (default estimator)"
```