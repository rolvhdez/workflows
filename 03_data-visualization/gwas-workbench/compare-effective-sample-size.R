# Load packages
suppressPackageStartupMessages(library(cli))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(R.utils))

# Load personalized functions
"%&%" <- function(a, b) paste0(a, b)
source("utils/00_utils.R")

# --------------------------------
# Read arguments
# --------------------------------

args <- commandArgs(trailingOnly = TRUE)
dir1 <- args[1]
dir2 <- args[2]
output_dir <- args[3]
model1 <- args[4]
model2 <- args[5]
model_label1 <- args[6]
model_label2 <- args[7]
pheno_dict <- args[8]

# Parameter error handling
if (!dir.exists(dir1)) cli_abort("{dir1} does not exist.")
if (!dir.exists(dir2)) cli_abort("{dir2} does not exist.")
if (!dir.exists(output_dir)) {
  # Provide a file in the directory for output
  output_dir <- output_dir %&% "."
}

# --------------------------------
# Read arguments
# --------------------------------
process_directory <- function(dir, model, label) {
  files <- list.files(dir, pattern = "*.gz", full.names = TRUE)
  df_list <- list()
  for (i in seq_along(files)) {
    f <- files[[i]]
    x <- process_sumstats(f, model, label)
    trait <- sub("^([0-9_]+).*", "\\1", basename(f))
    x$TRAIT <- trait
    df_list[[i]] <- x
  }
  df <- bind_rows(df_list)
  return(df)
}

df1 <- process_directory(dir1, model1, model_label1)
df2 <- process_directory(dir2, model2, model_label2)

#
traits <- intersect(unique(df1$TRAIT), unique(df2$TRAIT))
df1 <- subset(df1, TRAIT %in% traits)
df2 <- subset(df2, TRAIT %in% traits)
df_combined <- dplyr::bind_rows(df1, df2)

if (!is.null(pheno_dict)) {
  dict <- read.table(pheno_dict, sep = "\t", header = TRUE)

  df_combined <- df_combined %>%
    left_join(dict, by = "TRAIT") %>%
    select(-TRAIT) %>%
    rename(TRAIT = DESCRIPTION)
}

#
neff <- df_combined %>%
  group_by(DATASET, TRAIT) %>%
  summarise(
    median = median(N),
    ci_lower = quantile(N, 0.025), # CI 95% (alpha = 0.025)
    ci_upper = quantile(N, 0.975), # CI 95% (alpha = 0.025)
    .groups = "drop"
  )

p <- neff %>%
  ggplot() +
  geom_point(
    aes(
      x = TRAIT, y = median,
      color = DATASET, shape = DATASET
    ),
    size = 3,
    position = position_dodge(width = 0.25)
  ) +
  geom_errorbar(
    aes(
      x = TRAIT,
      ymin = ci_lower, ymax = ci_upper,
      color = DATASET
    ),
    width = 0.15,
    position = position_dodge(width = 0.25)
  ) +
  xlab("") +
  ylab(expression(Median~N[eff])) +
  labs(
    color = "Estimator",
    shape = "Estimator"
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(values = c(blue, red)) +
  theme(legend.position = "top") +
  coord_flip()

export_plot(p,
  output_dir %&% "effective-sample-sizes.png",
  width = 1920, height = 1920 * 0.60,
  res = 200, units = "px"
)