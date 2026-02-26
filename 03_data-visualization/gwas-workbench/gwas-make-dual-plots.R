# Author: Roberto Olvera Hernandez
# Date: 2025-11-07

# Load libraries
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

# --------------------------------
# Read arguments
# --------------------------------

args <- commandArgs(trailingOnly = TRUE)

sumstats1 <- args[1]
sumstats2 <- args[2]
output_dir <- args[3]
model1 <- if (length(args) < 4) "snipar" else args[4] # Default: snipar (v0.0.22)
model2 <- if (length(args) < 5) "snipar" else args[5] # Default: snipar (v0.0.22)
model_label1 <- if (length(args) < 6) "A" else args[6] # Default: snipar (v0.0.22)
model_label2 <- if (length(args) < 7) "B" else args[7] # Default: snipar (v0.0.22)
compute_bonferroni <- if (length(args) < 8) "no" else args[8] # Default: empty
phenotype <- if (length(args) < 9) "" else args[9] # Default: empty
annotations <- if (length(args) < 10) "no" else args[10] # Default: empty

# Parameter error handling
if (!file.exists(sumstats1)) cli_abort("{sumstats1} does not exist.")
if (!file.exists(sumstats2)) cli_abort("{sumstats2} does not exist.")
if (!dir.exists(output_dir)) {
  # Provide a file in the directory for output
  output_dir <- output_dir %&% "."
}

# --------------------------------
# Read arguments
# --------------------------------

df1 <- process_sumstats(sumstats1, model1, "a")
df2 <- process_sumstats(sumstats2, model2, "b")
df_combined <- dplyr::bind_rows(df1, df2)

if (compute_bonferroni == "yes") {
  bonferroni <- 0.05 / nrow(df_sumstats)
} else {
  bonferroni <- 5e-8 # Bonferroni adjusted P-Value
}

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

cli_alert_info(scales::comma(nrow(df1)) %&% " SNPs found in `" %&% sumstats1 %&% "`.")
cli_alert_info(scales::comma(nrow(df2)) %&% " SNPs found in `" %&% sumstats2 %&% "`.")
cli_alert_info("Bonferroni adjusted P-value: " %&% scales::scientific(bonferroni))

# --------------------------------
# QQ-Plots
# --------------------------------

qq1 <- get_qqvalues(df1)
qq2 <- get_qqvalues(df2)
pvalues <- list(qq1[[1]], qq2[[1]])
lambda <- list(qq1[[2]], qq2[[2]])

qqcaption <- paste0(
  "Lambda (", model_label1,"): ", lambda[[1]], "\n",
  "Lambda (", model_label2,"): ", lambda[[2]]
)

pvalues_combined <- rbind(
  cbind(pvalues[[1]], model = model_label1),
  cbind(pvalues[[2]], model = model_label2)
)

qqplot <- ggplot(data = pvalues_combined) +
  geom_point(
    aes(x = theoretical, y = observed, color = model),
    size = 1, shape = 16, alpha = 0.65
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_color_manual(values = c(
    model_label1 = blue,
    model_label2 = red
  )) +
  xlab(expression(Theoretical ~ -log[10](italic(p)))) +
  ylab(expression(Observed ~ -log[10](italic(p)))) +
  labs(
    title = phenotype,
    caption = qqcaption,
    color = "Model"
  )
export_plot(qqplot, output_dir %&% "dual-qqplot.png")
