simulateGZ <- function(N, fractions_cat1, fractions_cat2, target_cor, num_sim, plot_results = FALSE) {
  remaining_cat1 <- 1 - sum(fractions_cat1)
  remaining_cat2 <- 1 - sum(fractions_cat2)

  full_fractions_cat1 <- c(fractions_cat1, remaining_cat1)
  full_fractions_cat2 <- c(fractions_cat2, remaining_cat2)

  freq_array <- array(0, dim = c(num_sim, length(full_fractions_cat1), length(full_fractions_cat2)))

  for (i in 1:num_sim) {
    max_cor <- generateMaxCorrelatedData(N, fractions_cat1, fractions_cat2, negative = FALSE, rand = TRUE)

    result <- adjustCorrelation(max_cor$variable1, max_cor$variable2, target_cor)

    G <- result$variable1
    Z <- result$variable2

    freq_table <- table(factor(G, levels = 1:length(full_fractions_cat1)),
                        factor(Z, levels = 1:length(full_fractions_cat2))) / N
    freq_array[i, , ] <- as.matrix(freq_table)
  }

  mean_freqs <- apply(freq_array, c(2, 3), mean)

  plot_list <- NULL
  if (plot_results) {
    plot_list <- list()
    idx <- 1
    for (g in 1:length(full_fractions_cat1)) {
      for (z in 1:length(full_fractions_cat2)) {
        cum_means <- cumsum(freq_array[, g, z]) / seq_along(freq_array[, g, z])

        plot_data <- data.frame(Iteration = 1:num_sim, MeanFrequency = cum_means)

        plot_list[[idx]] <- ggplot(plot_data, aes(x = Iteration, y = MeanFrequency)) +
          geom_line(color = "blue") +
          labs(title = sprintf("G=%d, Z=%d", g, z),
               x = "Iteration", y = "Mean Frequency") +
          theme_minimal() +
          theme(plot.title = element_text(size = 8), axis.text = element_text(size = 6),
                axis.title = element_text(size = 7))
        idx <- idx + 1
      }
    }
  }

  list(MeanFrequencies = mean_freqs, FrequencyPlot = plot_list)
}



# Function to convert the mean frequencies matrix to a data frame
convertMeanFrequenciesToTable <- function(mean_freqs) {
  num_G <- dim(mean_freqs)[1]
  num_Z <- dim(mean_freqs)[2]

  freq_df <- data.frame(expand.grid(G = 1:num_G, Z = 1:num_Z),
                        q = as.vector(mean_freqs))

  freq_df$G <- freq_df$G
  freq_df$Z <- freq_df$Z

  colnames(freq_df) <- c("G", "Z", "q")

  return(freq_df)
}
