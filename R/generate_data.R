generate_data <- function(Beta0, Beta1, phi=c(),LD.r, P_g, P_z, N, n2, p_val, num_categories, model_type) {
  # set.seed(123)

  cat("Beta0 =", Beta0, " Beta1 =", Beta1, "\n")
  P_G <- 1 - P_g
  P_Z <- 1 - P_z
  LD <- LD.r * sqrt(P_G * P_g * P_Z * P_z)
  h.freqs <- rep(0, 4)
  h.freqs[1] <- P_G * P_Z + LD
  h.freqs[2] <- P_G * P_z - LD
  h.freqs[3] <- P_g * P_Z - LD
  h.freqs[4] <- P_g * P_z + LD

  if (any(h.freqs < 0)) {
    h.freqs <- pmax(h.freqs, 1e-05)
    h.freqs <- h.freqs / sum(h.freqs)
  }

  pval <- 1
  k <- 0
  while (pval >= p_val) {
    h1 <- sample(1:4, size = N, prob = h.freqs, replace = TRUE)
    h2 <- sample(1:4, size = N, prob = h.freqs, replace = TRUE)
    G <- as.numeric(h1 == 3 | h1 == 4) + as.numeric(h2 == 3 | h2 == 4)
    Z <- as.numeric(h1 == 2 | h1 == 4) + as.numeric(h2 == 2 | h2 == 4)

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
      for (j in 1:num_categories) {
        p_vals[, j] <- exp_sum_matrix[, j] / rowSums(exp_sum_matrix)
      }
      # p_vals[, 1] <- 1- rowSums(p_vals[, 2:(num_categories)])
    } else {
      stop("Invalid model type")
    }


    ylist <- vector("list", N)
    for (i in 1:N) {
      ydummy <- rmultinom(n = 1, size = 1, prob = p_vals[i, ])
      ylist[[i]] <- which(ydummy == 1)
    }

    Y <- unlist(ylist)


    if (model_type == "Proportional_Odds") {

      lmfit.comp <- vglm(Y ~ Z, family = cumulative(parallel = TRUE, reverse = FALSE))

    } else if (model_type == "Adjacent_Category") {

      lmfit.comp <- vglm(Y ~ Z, family = acat(reverse = TRUE, parallel = TRUE))

    } else if (model_type == "Stopping_Ratio") {

      lmfit.comp <- vglm(Y ~ Z, family = sratio(reverse = FALSE, parallel = TRUE))

    } else if (model_type == "Stereotype_Regression") {

      lmfit.comp <- rrvglm(-Y ~ Z, family = multinomial)

    }



    if (model_type == "Stereotype_Regression") {
      summary_fit <- VGAM:::summary.rrvglm(lmfit.comp)
    } else {
      summary_fit <- VGAM:::summaryvglm(lmfit.comp)
    }


    if (model_type == "Stereotype_Regression") {
      pval <- summary_fit@coef3[(2*num_categories-2),4]
    } else {
      pval <- summary_fit@coef3[num_categories, 4]
    }

    k <- k + 1
    if (pval < p_val) {
      G0 <- sample(0:2, size = N, prob = c(P_G^2, 2 * P_G * P_g, P_g^2), replace = TRUE)
      dat_sim_it <- data.frame(wait_it = k, Y = Y, G1 = G, Z = Z, G0 = G0)
    }
  }
  return(dat_sim_it)
}
