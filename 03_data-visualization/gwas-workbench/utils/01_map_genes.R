create_gene_ranges <- function() {
  # Download the Ensembl variant annotations ---
  ensembl_url <- "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz"
  ensembl_path <- "/tmp/Homo_sapiens.GRCh38.115.gtf.gz"
  if(!file.exists(ensembl_path)){
    cli_alert_info("Ensembl data base not found, downloading...")
    download.file(
      ensembl_url,
      destfile = ensembl_path,
      method = "wget",
      extra = "-r -p --random-wait",
      quiet = TRUE
    )
  }
  cli_alert_info("Reading " %&% ensembl_path %&% "...")
  df_ensembl <- as.data.frame(rtracklayer::import(ensembl_path))
  df_genes <- df_ensembl %>%
    filter(
      type == "gene",
      gene_biotype == "protein_coding"
    )
  gene_ranges <- GRanges(
    seqnames = df_genes$seqnames,
    IRanges(start = df_genes$start, end = df_genes$end)
  )
  names(gene_ranges) <- df_genes$gene_name
  return(gene_ranges)
}
annotate_genes_to_sig_snps <- function(sumstats, gene_range, sig = 5e-8){
  #' Function obtained from MCPS' GitHub organization
  #'
  #' @param sumstats Data frame (formated) with GWAS summary statistics.
  #' @param gene_range List of genetic ranges.
  #' @param sig Significance threshold. Default is genome-wide significance.

  require(dplyr)
  sumstats_sig <- dplyr::filter(sumstats, P <= sig)

  # Create the genetic ranges
  gene_nearest <- c()
  for (i in 1:dim(sumstats_sig)[1]) {
    granges_sig <- GRanges(
      seqnames = sumstats_sig$CHR[i],
      ranges = IRanges(start = sumstats_sig$BP[i], end = sumstats_sig$BP[i])
    )
    # Find the nearest gene
    g <- gene_range[(nearest(granges_sig, gene_range))] %>% names(.)
    gene_nearest <- append(gene_nearest, g)
  }
  sumstats_sig$GENE <- gene_nearest
  out_df <- c()
  for (gene in unique(sumstats_sig$GENE)){
    sub_df <- filter(sumstats_sig, GENE == gene) %>% arrange(P)
    out_df <- rbind(out_df, sub_df[1,])
  }
  return(out_df)
}
