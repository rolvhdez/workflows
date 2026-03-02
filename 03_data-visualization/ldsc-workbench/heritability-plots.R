#
# Author: Roberto Olvera Hernandez
# Date: 2026-02-12

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
h2_file <- args[1]
output_dir <- args[2]
pheno_file <- args[3]

# --------------------------------
df_h2 <- read.table(h2_file, header = TRUE, sep = "\t")
pheno_dict <- read.table(pheno_file, header = TRUE, sep = "\t")
df_h2 <- df_h2 %>%
  left_join(pheno_dict, by = "trait")

#
p <- df_h2 %>%
  ggplot(aes(x = reorder(description, h2), y = h2)) +
  geom_errorbar(aes(
    x = description,
    ymin = h2 - h2_se,
    ymax = h2 + h2_se
  ), width = 0.25) +
  geom_hline(yintercept = 0, color = red, linetype = "dashed") +
  geom_point(size = 3, shape = 16) +
  xlab("") +
  ylab(expression(h[SNP]^2)) +
  theme(
    legend.position = "none"
  ) +
  coord_flip()
export_plot(p,
  output_dir %&% "traits_heritability.png",
  height = 2 * 1080 * 0.75, width = 2 * 1080,
  res = 300, units = "px"
)