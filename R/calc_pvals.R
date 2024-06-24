calc_pvals <- function(Beta0, Beta1, model_type, num_categories, G, phi = NULL){
  if (model_type == "Proportional_Odds") {
    e_vals <- sapply(1:(num_categories - 1), function(x) exp(Beta0[x] + Beta1 * G))
    cum_probs <- e_vals / (1+e_vals)
    # Initialize the matrix for category-specific probabilities
    p_vals <- matrix(nrow = length(G), ncol = num_categories)
    p_vals[, 1] <- cum_probs[, 1]
    for (j in 2:(num_categories - 1)) {
      p_vals[, j] <- cum_probs[, j] - cum_probs[, j - 1]
    }
    p_vals[, num_categories] <- 1 - cum_probs[, num_categories - 1]
  } else if (model_type == "Adjacent_Category") {
    # Create a matrix with each row containing Beta0 values
    Beta0_matrix <- matrix(rep(Beta0, length(G)), nrow = length(G), ncol = (num_categories-1), byrow = TRUE)

    # Create a matrix with each row containing G values multiplied by Beta1
    G_matrix <- matrix(rep(G * Beta1,each=num_categories-1),  nrow = length(G), ncol=num_categories-1, byrow = T)

    # Sum the Beta0 matrix and G matrix
    sum_matrix <- Beta0_matrix + G_matrix

    # Calculate the cumulative sum for each row (from right to left)
    exp_cum_sum_matrix <- exp(t(apply(sum_matrix, 1, function(x) rev(cumsum(rev(x))))))


    p_vals <- matrix(nrow =length(G), ncol = num_categories)

    for (i in 1:length(G)) {
      # Compute denominator for the probabilities
      denominator <- 1 + sum(exp_cum_sum_matrix[i, ])

      # calculate probabilities for each category except the last one
      for (j in 1:(num_categories - 1)) {
        p_vals[i, j] <- exp_cum_sum_matrix[i, j] / denominator
      }

      # Calculate the probability for the last category
      p_vals[i, num_categories] <- 1 / denominator
    }
  } else if (model_type == "Stopping_Ratio") {
    Beta0_matrix <- matrix(rep(Beta0, length(G)), nrow = length(G), ncol = num_categories-1 , byrow = TRUE)

    # Create a matrix with each row containing G values multiplied by Beta1
    G_matrix <- matrix(rep(G * Beta1,each=num_categories-1),  nrow = length(G), ncol=num_categories-1, byrow = T)

    # Add the Beta0 matrix and G matrix
    exp_sum_matrix <- exp(Beta0_matrix + G_matrix)

    p_vals <- matrix(nrow = length(G), ncol = num_categories)

    # Iterate over each row to determine the probabilities
    for (i in 1:length(G)) {
      # Compute probabilities for each category
      for (j in 1:num_categories-1) {
        if (j == 1) {
          p_vals[i, j] <- exp_sum_matrix[i, j] / (1 + exp_sum_matrix[i, j])
        } else {
          product_term <- 1
          for (r in 1:(j-1)) {
            product_term <- product_term * (1 - exp_sum_matrix[i, r] / (1 + exp_sum_matrix[i, r]))
          }
          p_vals[i, j] <- exp_sum_matrix[i, j] / (1 + exp_sum_matrix[i, j]) * product_term
        }
      }
    }

    p_vals[, num_categories] <- 1 - rowSums(p_vals[, 1:(num_categories - 1)])

  } else if (model_type == "Stereotype_Regression") {
    Beta0_matrix <- matrix(rep(Beta0, length(G)), nrow = length(G), ncol = num_categories-1, byrow = TRUE)

    # Create a matrix with each row containing G values multiplied by Beta1

    G_matrix <- matrix(G * Beta1,ncol=1) %*% matrix(phi, nrow=1)
    Beta0_matrix <- matrix(rep(Beta0, length(G)), nrow = length(G), ncol = num_categories-1 , byrow = TRUE)
    # Add the Beta0 matrix and G matrix
    exp_sum_matrix <- cbind(rep(1,length(G)),exp(Beta0_matrix + G_matrix))
    # generate probabilities
    p_vals <- matrix(nrow = length(G), ncol = num_categories)
    for (j in 2:num_categories) {
      p_vals[, j] <- exp_sum_matrix[, j] / rowSums(exp_sum_matrix)
    }
    p_vals[, 1] <- 1- rowSums(p_vals[, 2:(num_categories)])
  }

  return(p_vals)
}
