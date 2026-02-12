make_effectsizes_plot <- function(df, title, bonferroni) {
  require(ggplot2)
  require(dplyr)

  k <- length(unique(df$SNP))
  caption <- paste0(
    "No. variants: ", scales::comma(k)
  )
  p <- ggplot() +
    # Significant variants
    geom_hline(
      aes(yintercept = -log10(bonferroni)), color = "red", linetype = "dashed"
    ) +
    geom_point(
      data = subset(df, df$P < bonferroni),
      aes(x = BETA, y = -log10(P), color = "Genome-wide significant (p < 5e-8)"),
    ) +
    geom_point(
      data = subset(df, df$P >= bonferroni),
      aes(x = BETA, y = -log10(P), color = "Other"),
      alpha = 0.85
    ) +
    # Non-Significant variants
    geom_hline(
      aes(yintercept = -log10(bonferroni)),
      color = red, linetype = "dashed"
    ) +
    xlab(expression(Effect~size~(beta))) +
    ylab(expression(-log[10](italic(p)))) +
    labs(
      title = title,
      caption = caption
    ) +
    scale_color_manual(
      name = "",
      values = c(
        "Other" = "gray",
        "Genome-wide significant (p < 5e-8)" = blue
      )
    ) +
    theme(
      legend.position = "top"
    )
}
make_effectmaf_plot <- function(df, title, bonferroni) {
  k <- length(unique(df$SNP))
  caption <- paste0(
    "No. variants: ", scales::comma(k)
  )
  p <- ggplot() +
    geom_point(data = subset(df, df$P > bonferroni),
      aes(x = MAF, y = BETA, color = "Other"),
      alpha = 0.85
    ) +
    geom_point(data = subset(df, df$P <= bonferroni),
      aes(x = MAF, y = BETA, color = "Genome-wide significant (p < 5e-8)"),
    ) +
    ylab(expression(Effect~size~(beta))) +
    xlab("Minor allele frequency (MAF)") +
    labs(
      title = title,
      caption = caption,
      color = expression(Effect~size~(beta))
    ) +
    scale_color_manual(
      name = "",
      values = c(
        "Other" = "gray",
        "Genome-wide significant (p < 5e-8)" = blue
      )
    ) +
    theme(
      legend.position = "top"
    )
  return(p)
}