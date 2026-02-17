# Author:
# Date: 
suppressPackageStartupMessages(library(cli))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

"%&%" <- function(a, b) paste0(a, b)
source("utils/00_utils.R")
source("utils/01_map_genes.R")

# Parameter ---
args <- commandArgs(trailingOnly = TRUE)
sumstats_a <- args[1]
sumstats_b <- args[2]
output_dir <- args[3]
phenotype <- if (length(args) < 4) "" else args[4] # Default: empty
model_a <- if (length(args) < 5) "snipar" else args[5]
model_b <- if (length(args) < 6) "snipar" else args[6]
model_label_a <- if (length(args) < 7) "Model A" else args[7]
model_label_b <- if (length(args) < 8) "Model B" else args[8]

# Check parameters ----
if (!file.exists(sumstats_a)) {
  print(sumstats_a)
  # Check if the summary statistics exist
  cli_abort(c(
    "{sumstats_a} does not exist",
    "x" = "You've supplied a file that does not exist."
  ))
}
if (!file.exists(sumstats_b)) {
  print(sumstats_b)
  # Check if the summary statistics exist
  cli_abort(c(
    "{sumstats_b} does not exist",
    "x" = "You've supplied a file that does not exist."
  ))
}
if (!dir.exists(output_dir)) {
  # Provide a file in the directory for output
  output_dir <- output_dir %&% "."
}

#------------------------------------

# Read the summary statistics -------
df_a <- fancy_process(
  process = read_sumstats_file,
  message = "Reading " %&% sumstats_a,
  # Function parameters
  sumstats_path = sumstats_a,
  chunk_size = 1000000
)
df_a <- reformat_sumstats(df_a, model_a)
df_a$MODEL <- "model_a"

#
df_b <- fancy_process(
  process = read_sumstats_file,
  message = "Reading " %&% sumstats_b,
  # Function parameters
  sumstats_path = sumstats_b,
  chunk_size = 1000000
)
df_b <- reformat_sumstats(df_b, model_b)
df_b$MODEL <- "model_b"

#
df_sumstats <- dplyr::bind_rows(df_a, df_b)

# #
# bonferroni <- 5e-8
# # Check that there are significant SNPs to annotate
# sig_k <- df_sumstats %>% filter(P <= bonferroni) %>% pull(SNP)
# if (length(sig_k) > 0) {
#   cli::cli_alert_warning(scales::comma(length(sig_k)) %&% " SNPs found at p <= " %&% bonferroni)
#   genes <- create_gene_ranges()
#   df_annotations <- annotate_genes_to_sig_snps(
#     sumstats = df_sumstats,
#     gene_range = genes,
#     sig = bonferroni
#   )
#   df_annotations <- df_annotations %>% select(SNP, GENE)
#   df_sumstats <- df_sumstats %>%
#     left_join(df_annotations, by = "SNP")
# } else {
#   cli::cli_alert_warning("No significant SNPs were found at p <= " %&% bonferroni %&% ". Skipping annotation.")
# }

# Correlations ---
effects <- df_sumstats %>%
  tidyr::pivot_wider(
    id_cols = SNP,
    names_from = MODEL,
    values_from = BETA
  ) %>%
  tidyr::drop_na()
pvalues <- df_sumstats %>%
  tidyr::pivot_wider(
    id_cols = SNP,
    names_from = MODEL,
    values_from = P
  ) %>%
  tidyr::drop_na()
effect_pval <- effects %>%
  left_join(pvalues, by = "SNP") %>%
  rename_with(~ gsub("\\.x$", "_effect", .)) %>%
  rename_with(~ gsub("\\.y$", "_pval", .))

r_value <- cor(effects$model_a, effects$model_b, method = "pearson")
r_value_msg <- paste0(
  "Pearson's r^2 = ", round(r_value, 4)
)
cli::cli_alert_info(r_value_msg)

#
# Create conditional axis labels
x_label <- if(tolower(model_a) == "snipar") {
  bquote(delta[.(model_label_a)])
} else {
  bquote(beta[.(model_label_a)])
}
y_label <- if(tolower(model_b) == "snipar") {
  bquote(delta[.(model_label_b)])
} else {
  bquote(beta[.(model_label_b)])
}

p <- effect_pval %>%
  mutate(
    significance = case_when(
      model_a_pval < 5e-8 & model_b_pval >= 5e-8 ~ model_label_a,
      model_b_pval < 5e-8 & model_a_pval >= 5e-8 ~ model_label_b,
      model_a_pval < 5e-8 & model_b_pval < 5e-8 ~ "Both",
      TRUE ~ "Neither"
    )
  ) %>%
  ggplot(aes(x = model_a_effect, y = model_b_effect)) +
  # First layer: All "Neither" points with transparency
  geom_point(
    data = . %>% filter(significance == "Neither"),
    aes(fill = significance),
    shape = 21,
    color = "gray",
    size = 1,
    alpha = 0.5
  ) +
  # Second layer: Significant points (no transparency)
  geom_point(
    data = . %>% filter(significance != "Neither"),
    aes(fill = significance),
    shape = 21,
    color = "black",
    size = 2
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text",
    x = 0,
    y = max(effect_pval$model_b_effect),
    label = paste("r = ", round(r_value, 4)),
    hjust = 0, vjust = 1, size = 3, color = "black"
  ) +
  scale_fill_manual(
    values = c(
      setNames(blue, model_label_a),
      setNames(green, model_label_b),
      "Both" = red,
      "Neither" = "lightgray"
    ),
    name = "Is significant (p < 5e-8)"
  ) +
  labs(title = phenotype) +
  xlab(x_label) +
  ylab(y_label) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 10)
  )
export_plot(p,
  output_dir %&% "effect_scatter." %&% model_a %&% "-" %&% model_b %&% ".png",
  height = 1080 * 0.75 * 2, width = 1080 * 2,
  res = 300, units = "px"
)
