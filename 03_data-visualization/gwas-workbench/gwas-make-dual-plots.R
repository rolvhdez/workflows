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
suppressPackageStartupMessages(library(patchwork))

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
phenotype <- if (length(args) < 8) "" else args[8] # Default: empty
annotations <- if (length(args) < 9) "no" else args[9] # Default: empty

# Parameter error handling
if (!file.exists(sumstats1)) cli_abort("{sumstats1} does not exist.")
if (!file.exists(sumstats2)) cli_abort("{sumstats2} does not exist.")
if (!dir.exists(output_dir)) {
  # Provide a file in the directory for output
  output_dir <- output_dir %&% "."
}

get_annotations <- function(df, bonferroni) {
  snps <- subset(df, P <= bonferroni)$SNP
  if (length(snps) > 0) {
    cli::cli_alert_warning(
      paste0(scales::comma(length(snps)), " SNPs found at p <= ", bonferroni)
    )
    genes <- create_gene_ranges()
    annotations <- annotate_genes_to_sig_snps(
      sumstats = df,
      gene_range = genes,
      sig = bonferroni
    )
    return(dplyr::select(annotations, SNP, GENE))
  } else {
    cli::cli_alert_warning(
      paste0("No significant SNPs were found at p <= ", bonferroni, ". Skipping annotation.")
    )
    return(NULL)
  }
}

# --------------------------------
# Read arguments
# --------------------------------

df1 <- process_sumstats(sumstats1, model1, "a")
df2 <- process_sumstats(sumstats2, model2, "b")
df_combined <- dplyr::bind_rows(df1, df2)

bonferroni <- 5e-8 # Bonferroni adjusted P-Value
cli_alert_info("Bonferroni adjusted P-value: " %&% scales::scientific(bonferroni))

if (annotations == "yes") {
  df1_annotations <- get_annotations(df1, bonferroni)
  df2_annotations <- get_annotations(df2, bonferroni)
  #
  if (!is.null(df1_annotations)) df1 <- df1 %>% left_join(df1_annotations, by = "SNP")
  if (!is.null(df2_annotations)) df2 <- df2 %>% left_join(df2_annotations, by = "SNP")
}

# --------------------------------
# QQ-Plots
# --------------------------------
qq1 <- get_qqvalues(df1)
qq2 <- get_qqvalues(df2)
pvalues <- list(qq1[[1]], qq2[[1]])
lambda <- list(qq1[[2]], qq2[[2]])

qqcaption <- paste0(
  "Lambda (", model_label1, "): ", round(lambda[[1]], 4), "\n",
  "Lambda (", model_label2, "): ", round(lambda[[2]], 4)
)
pvalues_combined <- rbind(
  cbind(pvalues[[1]], model = model_label1),
  cbind(pvalues[[2]], model = model_label2)
)
qqplot <- ggplot(data = pvalues_combined) +
  geom_point(
    aes(x = theoretical, y = observed, color = model),
    size = 2, shape = 16, alpha = 0.65
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_color_manual(values = setNames(
    c(blue, red),
    c(model_label1, model_label2)
  )) +
  xlab(expression(Theoretical ~ -log[10](italic(p)))) +
  ylab(expression(Observed ~ -log[10](italic(p)))) +
  labs(
    title = phenotype,
    caption = qqcaption,
    color = "Model"
  ) +
  theme(legend.position = "top")
export_plot(qqplot, output_dir %&% "dual-qqplot.png")

# --------------------------------
# Manhattan Plot
# --------------------------------
process_manhattan <- function(df, bonferroni, phenotype) {
  m_list <- format_for_manhattan(df)
  df_m <- m_list[[1]]
  df_a <- m_list[[2]]
  m <- make_manhattan(df_m, df_a, "", bonferroni)
  return(m)
}

m1 <- process_manhattan(df1, bonferroni, phenotype) +
  labs(x = "") +
  ggtitle(model_label1) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 8))
if ("GENE" %in% names(df1)) {
  m1 <- m1 +
    geom_text_repel(
      data = subset(df1, !is.na(df1$GENE)),
      aes(label = GENE),
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = 16,
      size = 3
    )
}

m2 <- process_manhattan(df2, bonferroni, phenotype) +
  ggtitle(model_label2) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 8))
if ("GENE" %in% names(df2)) {
  m2 <- m2 +
    geom_text_repel(
      data = subset(df2, !is.na(df2$GENE)),
      aes(label = GENE),
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = 16,
      size = 3
    )
}

manhattan_plot <- m1 / m2 +
  plot_annotation(title = phenotype)
export_plot(manhattan_plot, output_dir %&% "dual-manhattan.png")