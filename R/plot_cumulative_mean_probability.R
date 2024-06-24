plot_cumulative_mean_probability <- function(data, category, G_value, true_probs_mat) {

  prob_column <- paste0("p", category)


  cumulative_means <- data %>%
    dplyr::filter(G == G_value) %>%
    dplyr::group_by(model_type) %>%
    dplyr::mutate(CumulativeMean = cumsum(get(prob_column)) / seq_along(get(prob_column))) %>%
    dplyr::ungroup() %>%
    dplyr::select(iteration, model_type, CumulativeMean)

  true_prob <- true_probs_mat[G_value, category]

  ggplot(cumulative_means, aes(x = iteration, y = CumulativeMean, color = model_type, group = model_type)) +
    geom_line() +
    geom_hline(yintercept = true_prob, linetype = "dashed", color = "black", size = 1) +
    labs(title = paste("Cumulative Mean of Probability for Category", category, "and G =", G_value),
         x = "Iteration", y = "Cumulative Mean Probability") +
    theme_minimal() +
    scale_color_brewer(palette = "Set1") +
    guides(color = guide_legend(title = "Model Type"))
}
