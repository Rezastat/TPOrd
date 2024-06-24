calc_stats <- function(matrix) {
  apply(matrix, 2, function(x) {
    c(mean = mean(x, na.rm = TRUE), variance = var(x, na.rm = TRUE), quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
  })
}
