format_for_manhattan <- function(sumstats){
  x <- sumstats %>%
    group_by(CHR) %>%
    summarise(CHR_LEN = max(BP)) %>%
    mutate(COORD = cumsum(as.numeric(CHR_LEN)) - CHR_LEN) %>%
    dplyr::select(-CHR_LEN) %>%
    right_join(sumstats, ., by = "CHR") %>%
    arrange(CHR, BP) %>%
    mutate(BP_CUM = COORD + BP)
  y <- x %>%
    dplyr::select(CHR, BP_CUM) %>%
    group_by(CHR) %>%
    summarise(CENTER = min(BP_CUM) + (max(BP_CUM) - min(BP_CUM)) / 2)
  return(list(manhattan = x, axis = y))
}
make_manhattan <- function(df, axis, title, bonferroni) {
  k <- length(unique(df$SNP))
  caption <- paste0(
    "No. variants: ", scales::comma(k), "\n",
    "Bonferroni adjusted p-value < " %&% scales::scientific(bonferroni, digits = 4) %&% " (red line)"
  )
  m <- df %>%
    ggplot(aes(x = BP_CUM, y = -log10(P))) +
    geom_hline( # High
      yintercept = -log10(bonferroni),
      color = red,
      linetype = "dashed",
      alpha = 0.8
    ) +
    geom_point(
      aes(color = as.factor(CHR)),
      shape = 16, alpha = 0.85
    ) +
    # Axis labels
    xlab("Chromosome") +
    ylab(expression(-log[10](italic(p)))) +
    labs(
      title = title,
      caption = caption
    ) +
    # Modify axis
    scale_x_continuous(
      label = df_axis$CHR,
      breaks = df_axis$CENTER
    ) +
    scale_color_manual(
      values = rep(c("gray", "steelblue"), 22)
    ) +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(face = "bold"),
      panel.spacing = unit(1.5, "lines"),
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5, hjust = 1
      )
    )
  return(m)
}