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
output <- args[2]

# --------------------------------
# 2.
# --------------------------------

df_h2 <- read.table(h2_file, header = TRUE, sep = "\t")
print(head(df_h2))