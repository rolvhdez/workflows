get_qqvalues <- function(sumstats, df = 1) {
  #' Compute the expected values of a linear model
  #' 
  #' @param sumstats Data frame with summary statistics
  p_values <- sumstats$P
  chi <- qchisq(1 - p_values, df = df)

  # Lambda inflation correction
  lambda <- median(chi, na.rm = TRUE) / qchisq(0.5, df)
  p_values <- chi / lambda
  n <- length(p_values)
  probs <- (1:n - 0.5) / n
  theoretical <- qchisq(probs, df)
  df <- data.frame(
    observed = sort(p_values),
    theoretical = sort(theoretical)
  )
  return(list(qqvalues = df, lambda = lambda))
}
make_qqplot <- function(pvalues, title, lambda){
  #' Make a QQ Plot for observed vs. theoretical p-values
  #' 
  #' @param pvalues A data frame with the theoretical and observed p-values
  #' @param title Title for the plot
   
  k <- nrow(pvalues)

  caption <- paste0(
    "No. variants: ", scales::comma(k), "\n",
    "Lambda inflation factor: ", round(lambda, 4)
  )

  x <- pvalues %>%
    filter(observed != Inf) %>%
    ggplot(aes(x = theoretical, y = observed)) +
    geom_point(size = 2, shape = 16, alpha = 0.65, color = "black") +
    geom_abline(
      slope = 1 * lambda, intercept = 0,
      color = red, linewidth = 1, alpha = 0.85,
      linetype = "dashed"
    ) +
    xlab(expression(Theoretical ~ -log[10](italic(p)))) +
    ylab(expression(Observed ~ -log[10](italic(p)))) +
    labs(
      title = title,
      caption = caption
    )
  return(x)
}