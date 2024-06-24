
sim_categorical_data <- function(Beta0, Beta1, phi=c(),target_cor, fractions_cat1, fractions_cat2, N, n2, cor_YZ_break, num_categories, model_type, negative = F, rand = T) {


  cat("Beta0 =", Beta0, " Beta1 =", Beta1, "\n")
  if (target_cor < 0) {
    negative = TRUE
  } else {
    negative = FALSE
  }
  best_cor_Y_Z <- 0
  best_result <- NULL
  condition = TRUE
  k <- 0
  while (condition & k<=200) {

    k <- k + 1

    max_cor <- generateMaxCorrelatedData(N, fractions_cat1, fractions_cat2, negative = negative, rand = rand)
    result <- adjustCorrelation(max_cor$variable1, max_cor$variable2, target_cor)
    G <- result$variable1
    Z <- result$variable2

    if (model_type == "Proportional_Odds") {
      e_vals <- sapply(1:(num_categories - 1), function(x) exp(Beta0[x] + Beta1 * G))
      cum_probs <- e_vals / (1+e_vals)
      # Initialize the matrix for category-specific probabilities
      p_vals <- matrix(nrow = N, ncol = num_categories)
      p_vals[, 1] <- cum_probs[, 1]
      for (j in 2:(num_categories - 1)) {
        p_vals[, j] <- cum_probs[, j] - cum_probs[, j - 1]
      }
      p_vals[, num_categories] <- 1 - cum_probs[, num_categories - 1]
    } else if (model_type == "Adjacent_Category") {
      # Create a matrix with each row containing Beta0 values
      Beta0_matrix <- matrix(rep(Beta0, N), nrow = N, ncol = (num_categories-1), byrow = TRUE)

      # Create a matrix with each row containing G values multiplied by Beta1
      G_matrix <- matrix(rep(G * Beta1,each=num_categories-1),  nrow = N, ncol=num_categories-1, byrow = T)

      # Sum the Beta0 matrix and G matrix
      sum_matrix <- Beta0_matrix + G_matrix

      # Calculate the cumulative sum for each row (from right to left)
      exp_cum_sum_matrix <- exp(t(apply(sum_matrix, 1, function(x) rev(cumsum(rev(x))))))


      p_vals <- matrix(nrow = N, ncol = num_categories)

      for (i in 1:N) {
        # Compute denominator for the probabilities
        denominator <- 1 + sum(exp_cum_sum_matrix[i, ])

        # Compute probabilities for each category except the last one
        for (j in 1:(num_categories - 1)) {
          p_vals[i, j] <- exp_cum_sum_matrix[i, j] / denominator
        }

        # Compute probability for the last category
        p_vals[i, num_categories] <- 1 / denominator
      }
    } else if (model_type == "Stopping_Ratio") {
      Beta0_matrix <- matrix(rep(Beta0, N), nrow = N, ncol = num_categories-1 , byrow = TRUE)

      # Create a matrix with each row containing G values multiplied by Beta1
      G_matrix <- matrix(rep(G * Beta1,each=num_categories-1),  nrow = N, ncol=num_categories-1, byrow = T)

      # Sum the Beta0 matrix and G matrix
      exp_sum_matrix <- exp(Beta0_matrix + G_matrix)

      p_vals <- matrix(nrow = N, ncol = num_categories)

      # Iterate over each row to compute probabilities
      for (i in 1:N) {
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
      Beta0_matrix <- matrix(rep(Beta0, N), nrow = N, ncol = num_categories-1, byrow = TRUE)

      # Create a matrix with each row containing G values multiplied by Beta1

      G_matrix <- matrix(G * Beta1,ncol=1) %*% matrix(phi, nrow=1)
      Beta0_matrix <- matrix(rep(Beta0, N), nrow = N, ncol = num_categories-1 , byrow = TRUE)
      # Sum the Beta0 matrix and G matrix
      exp_sum_matrix <- cbind(rep(1,N),exp(Beta0_matrix + G_matrix))
      # generate probabilities
      p_vals <- matrix(nrow = N, ncol = num_categories)
      for (j in 2:num_categories) {
        p_vals[, j] <- exp_sum_matrix[, j] / rowSums(exp_sum_matrix)
      }
      p_vals[, 1] <- 1- rowSums(p_vals[, 2:(num_categories)])
    } else {
      stop("Invalid model type")
    }


    ylist <- vector("list", N)
    for (i in 1:N) {
      ydummy <- rmultinom(n = 1, size = 1, prob = p_vals[i, ])
      ylist[[i]] <- which(ydummy == 1)
    }

    Y <- unlist(ylist)

    cor_Y_G <- abs(cor(Y, G, method = "spearman"))
    cor_Y_Z <- abs(cor(Y, Z, method = "spearman"))

    if (cor_Y_Z > best_cor_Y_Z) {
      if (cor_Y_Z > cor_Y_G * cor_YZ_break) {
        condition <- F
      }
      best_cor_Y_Z <- cor_Y_Z
      G0 <- sample(1:max(unique(G)), size = N, prob = table(G)/N, replace = TRUE) # Generating G0 based on G
      # Update best_result with current iteration's data
      best_result <- list(Y = Y, G1 = G, Z = Z, G0 = G0, cor_Y_Z = cor_Y_Z, cor_Y_G = cor_Y_G, wait_it = k)


    }
  }

  # Prepare the output based on the best_result
  if (!is.null(best_result)) {
    dat_sim_it <- data.frame(wait_it = best_result$wait_it, Y = best_result$Y, G1 = best_result$G1, Z = best_result$Z, G0 = best_result$G0)
  } else {
    dat_sim_it <- data.frame(wait_it = integer(), Y = integer(), G1 = numeric(), Z = numeric(), G0 = numeric())
    warning("No iterations met the criteria, returning an empty data frame.")
  }

  return(dat_sim_it)
}




#   function(Beta0, Beta1, phi=c(),target_cor, fractions_cat1, fractions_cat2, N, n2, p_val, num_categories, model_type, negative = F, rand = T) {
#
#
#   cat("Beta0 =", Beta0, " Beta1 =", Beta1, "\n")
#   if (target_cor < 0) {
#     negative = TRUE
#   } else {
#     negative = FALSE
#   }
# # initiate p
#   pval <- 1
#   k <- 0
#   while (pval >= p_val) {
#     max_cor <- generateMaxCorrelatedData(N, fractions_cat1, fractions_cat2, negative = negative, rand = rand)
#     result <- adjustCorrelation(max_cor$variable1, max_cor$variable2, target_cor)
#     G <- result$variable1
#     Z <- result$variable2
#
#     if (model_type == "Proportional_Odds") {
#       e_vals <- sapply(1:(num_categories - 1), function(x) exp(Beta0[x] + Beta1 * G))
#       cum_probs <- e_vals / (1+e_vals)
#       # Initialize the matrix for category-specific probabilities
#       p_vals <- matrix(nrow = N, ncol = num_categories)
#       p_vals[, 1] <- cum_probs[, 1]
#       for (j in 2:(num_categories - 1)) {
#         p_vals[, j] <- cum_probs[, j] - cum_probs[, j - 1]
#       }
#       p_vals[, num_categories] <- 1 - cum_probs[, num_categories - 1]
#     } else if (model_type == "Adjacent_Category") {
#       # Create a matrix with each row containing Beta0 values
#       Beta0_matrix <- matrix(rep(Beta0, N), nrow = N, ncol = (num_categories-1), byrow = TRUE)
#
#       # Create a matrix with each row containing G values multiplied by Beta1
#       G_matrix <- matrix(rep(G * Beta1,each=num_categories-1),  nrow = N, ncol=num_categories-1, byrow = T)
#
#       # Sum the Beta0 matrix and G matrix
#       sum_matrix <- Beta0_matrix + G_matrix
#
#       # Calculate the cumulative sum for each row (from right to left)
#       exp_cum_sum_matrix <- exp(t(apply(sum_matrix, 1, function(x) rev(cumsum(rev(x))))))
#
#
#       p_vals <- matrix(nrow = N, ncol = num_categories)
#
#       for (i in 1:N) {
#         # Compute denominator for the probabilities
#         denominator <- 1 + sum(exp_cum_sum_matrix[i, ])
#
#         # Compute probabilities for each category except the last one
#         for (j in 1:(num_categories - 1)) {
#           p_vals[i, j] <- exp_cum_sum_matrix[i, j] / denominator
#         }
#
#         # Compute probability for the last category
#         p_vals[i, num_categories] <- 1 / denominator
#       }
#     } else if (model_type == "Stopping_Ratio") {
#       Beta0_matrix <- matrix(rep(Beta0, N), nrow = N, ncol = num_categories-1 , byrow = TRUE)
#
#       # Create a matrix with each row containing G values multiplied by Beta1
#       G_matrix <- matrix(rep(G * Beta1,each=num_categories-1),  nrow = N, ncol=num_categories-1, byrow = T)
#
#       # Sum the Beta0 matrix and G matrix
#       exp_sum_matrix <- exp(Beta0_matrix + G_matrix)
#
#       p_vals <- matrix(nrow = N, ncol = num_categories)
#
#       # Iterate over each row to compute probabilities
#       for (i in 1:N) {
#         # Compute probabilities for each category
#         for (j in 1:num_categories-1) {
#           if (j == 1) {
#             p_vals[i, j] <- exp_sum_matrix[i, j] / (1 + exp_sum_matrix[i, j])
#           } else {
#             product_term <- 1
#             for (r in 1:(j-1)) {
#               product_term <- product_term * (1 - exp_sum_matrix[i, r] / (1 + exp_sum_matrix[i, r]))
#             }
#             p_vals[i, j] <- exp_sum_matrix[i, j] / (1 + exp_sum_matrix[i, j]) * product_term
#           }
#         }
#       }
#
#       p_vals[, num_categories] <- 1 - rowSums(p_vals[, 1:(num_categories - 1)])
#
#     } else if (model_type == "Stereotype_Regression") {
#       Beta0_matrix <- matrix(rep(Beta0, N), nrow = N, ncol = num_categories-1, byrow = TRUE)
#
#       # Create a matrix with each row containing G values multiplied by Beta1
#
#       G_matrix <- matrix(G * Beta1,ncol=1) %*% matrix(phi, nrow=1)
#       Beta0_matrix <- matrix(rep(Beta0, N), nrow = N, ncol = num_categories-1 , byrow = TRUE)
#       # Sum the Beta0 matrix and G matrix
#       exp_sum_matrix <- cbind(rep(1,N),exp(Beta0_matrix + G_matrix))
#       # generate probabilities
#       p_vals <- matrix(nrow = N, ncol = num_categories)
#       for (j in 2:num_categories) {
#         p_vals[, j] <- exp_sum_matrix[, j] / rowSums(exp_sum_matrix)
#       }
#       p_vals[, 1] <- 1- rowSums(p_vals[, 2:(num_categories)])
#     } else {
#       stop("Invalid model type")
#     }
#
#
#     ylist <- vector("list", N)
#     for (i in 1:N) {
#       ydummy <- rmultinom(n = 1, size = 1, prob = p_vals[i, ])
#       ylist[[i]] <- which(ydummy == 1)
#     }
#
#     Y <- unlist(ylist)
#
#
#     fitResult <- capture_fit_with_warnings(Y = Y, Z = Z, model_type = model_type)
#
#     if (is.null(fitResult$lmfit.comp) || fitResult$warningOccurred) {
#       next
#     } else {
#       lmfit.comp <- fitResult$lmfit.comp
#
#     if (model_type == "Stereotype_Regression") {
#       wald_statistic <- (lmfit.comp[["beta"]][1,1])^2 / (lmfit.comp[["beta"]][1,2])^2
#     } else {
#       summary_fit <- VGAM:::summaryvglm(lmfit.comp)
#     }
#
#
#     if (model_type == "Stereotype_Regression") {
#       pval <- pchisq(wald_statistic, df = 1, lower.tail = FALSE)
#     } else {
#       pval <- summary_fit@coef3[num_categories, 4]
#     }
#
#     k <- k + 1
#     if (pval < p_val) {
#       G0 <- sample(1:max(unique(result$variable1)), size = N, prob = table(result$variable1)/N, replace = TRUE)
#       dat_sim_it <- data.frame(wait_it = k, Y = Y, G1 = G, Z = Z, G0 = G0)
#       break # Assuming you stop once a condition is met; remove if you want to continue the loop
#     }
#     }
#   }
#
#   return(dat_sim_it)
# }
#
#
#
#
#
#
#
#
#
