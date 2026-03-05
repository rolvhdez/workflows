rm(list=ls())
library(dplyr)
library(scales)
library(cli)

# 
'%&%' <- function(a, b) paste0(a, b)
get_firstdegree_ids <- function(df) {
  require(dplyr)
  mask <- df$InfType %in% c('PO', 'FS', 'Dup/MZ')
  df_new <- df[mask, ]
  ids <- c(df_new$ID1, df_new$ID2)
  ids <- unique(ids)
  return(ids)
}

set.seed(123)

# Obtain the individuals with first degree relatives ----
king <- read.table("~/data/king_ibdseg_4th.seg", sep = "\t", head = TRUE)
baseline <- read.table("~/data/complete-set.qc-phenotypes.txt", sep = "\t", head = TRUE)

# Create phenotype filter ---
first_degree_ids <- get_firstdegree_ids(king)
extra_firstdegree <- sample(first_degree_ids, ceiling(length(first_degree_ids) / 3))
baseline_mask <- (
  baseline$IID %in% extra_firstdegree | !baseline$IID %in% first_degree_ids
) & !grepl("-DUP$", baseline$IID)
unrelated_baseline <- baseline[baseline_mask, ]
head(unrelated_baseline)

# Remove the 'X' from the headers ---
pheno_colnames <- names(unrelated_baseline)[15:length(unrelated_baseline)]
pheno_colnames <- gsub("^X", "", pheno_colnames)
names(unrelated_baseline)[15:length(unrelated_baseline)] <- pheno_colnames

# 
num_ids <- length(unique(baseline$IID))
cli_alert(
  paste0(
  'Individuals with first-degree relatives: ',
  scales::comma(length(first_degree_ids)),
  ' (', scales::percent(length(first_degree_ids) / num_ids, 2), ')'
  )
)
unrelated_ids_num <- length(unique(unrelated_baseline$IID))
cli_alert(paste0(
  "Individuals remaining: ",
  scales::comma(unrelated_ids_num),
  " (", scales::percent( unrelated_ids_num / num_ids ), ")"
))
write.table(
  unrelated_baseline,
  '~/data/complete-set.unrelated.txt',
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
