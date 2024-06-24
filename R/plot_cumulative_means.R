plot_cumulative_means <- function(coefs_matrix, title = "Cumulative Mean of Estimates") {

  cumulative_means <- apply(coefs_matrix, 2, function(x) cumsum(x) / seq_along(x))


  cumulative_means_df <- as.data.frame(cumulative_means)
  cumulative_means_long <- reshape2::melt(cumulative_means_df, variable.name = "Coefficient", value.name = "CumulativeMean")
  cumulative_means_long$Iteration <- rep(1:nrow(cumulative_means_df), ncol(cumulative_means_df))


  pl <- ggplot(cumulative_means_long, aes(x = Iteration, y = CumulativeMean, color = Coefficient)) +
    geom_line() +
    labs(title = title, x = "Iteration", y = "Cumulative Mean") +
    theme_minimal() +
    scale_color_discrete(name = "Coefficient")

  print(pl)
}
