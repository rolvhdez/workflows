#
# Author: Roberto Olvera Hernandez
# Date: 2026-02-12

#
# --------------------------------
# 1. Set up
# --------------------------------
# Load libraries
suppressPackageStartupMessages(library(cli))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

# Load personalized functions
"%&%" <- function(a, b) paste0(a, b)
source("utils/00_utils.R")

# Process arguments
args <- commandArgs(trailingOnly = TRUE)
h2_a_file <- args[1]
h2_b_file <- args[2]
output_dir <- args[3]
pheno_file <- args[4]
label_a <- args[5]
label_b <- args[6]

#
df_h2_a <- read.table(h2_a_file, header = TRUE, sep = "\t")
df_h2_a$model <- "a"

df_h2_b <- read.table(h2_b_file, header = TRUE, sep = "\t")
df_h2_b$model <- "b"
df_h2 <- dplyr::bind_rows(df_h2_a, df_h2_b)

pheno_dict <- read.table(pheno_file, header = TRUE, sep = "\t")

#
df_long <- df_h2 %>%
  tidyr::pivot_wider(
    id_cols = trait,
    names_from = model,
    values_from = c(h2, h2_se)
  ) %>%
  tidyr::drop_na() %>%
  left_join(pheno_dict, by = "trait") %>%
  mutate(ratio = h2_b / h2_a)

#
x_label <- bquote(h[SNP]^2~(.(label_a)))
y_label <- bquote(h[SNP]^2~(.(label_b)))
p <- df_long %>%
  ggplot(aes(x = h2_a, y = h2_b)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "lightgray"
  ) +
  geom_errorbar(aes(
    x = h2_a,
    ymin = h2_b - (h2_se_b * 1.96),
    ymax = h2_b + (h2_se_b * 1.96)
  ), color = "lightgray", width = 0.005) +
  geom_errorbarh(aes(
    y = h2_b,
    xmin = h2_a - (h2_se_a * 1.96),
    xmax = h2_a + (h2_se_a * 1.96)
  ), color = "lightgray", width = 0.005) +
  geom_point(
    aes(fill = ratio),
    size = 2.55, shape = 21
  ) +
  scale_fill_gradient2(
    low = blue, high = red,
    mid = "white", midpoint = 1
  ) +
  ggrepel::geom_text_repel(aes(
    x = h2_a, y = h2_b,
    label = description
  ),
    size = 3, box.padding = 0.65,
    point.padding = 0.65,
    max.overlaps = 16
  ) +
  xlab(x_label) +
  ylab(y_label) +
  theme(
    legend.position = "none"
  )

export_plot(p,
  output_dir %&% "heritability_comparison.png",
  height = 2 * 1080 * 0.75, width = 2 * 1080,
  res = 300, units = "px"
)