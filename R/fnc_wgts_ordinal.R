fnc_wgts_ordinal <- function(formula, family, dat0, theta, q, G_, Z_, uniqZ, num_categories,N_total,n2, model_type) {
  wg0 <- id_ <- NULL

  resp <- as.character(formula[[2]])
  terms_q <- c(G_, Z_)
  pG_ <- aggregate(as.formula(paste0("q~", paste(G_, collapse = "+"))), q, FUN = sum)
  names(pG_)[ncol(pG_)] <- "q_g"

  dat0$ord <- 1:nrow(dat0)
  dat0bis <- merge(as.data.table(dat0), as.data.table(q[, c(terms_q, "q")]), by = terms_q, all.x = TRUE, sort = FALSE)
  dat0bis <- merge(dat0bis, as.data.table(pG_), by = G_, all.x = TRUE, sort = FALSE)

  dat0bis$ind <- ifelse(dat0bis[[Z_]] %in% uniqZ, 1, 0)
  dat0bis$q <- dat0bis$ind * dat0bis$q + (1 - dat0bis$ind) * dat0bis$q_g

  X0 = model.matrix(formula, dat0bis)
  extended_X0 <-  expand_model_matrix(X0, num_categories-1)


  n2_bar <-(N_total-n2) * max(unique(q[,"G"]))


  if (model_type == "Stereotype_Regression"){
    phi_col <- rep(c(theta[1:(num_categories-2)],1),n2_bar)
    extended_X0[, num_categories:ncol(extended_X0)] <- extended_X0[, num_categories:ncol(extended_X0)] * phi_col
    eta_vals <- extended_X0 %*% theta[(num_categories-1):length(theta)]
  } else {
    eta_vals <- extended_X0 %*% theta
  }
  eta_matrix <- matrix(eta_vals,nrow= n2_bar,ncol = num_categories-1, byrow = T)

  if (model_type == "Proportional_Odds") {

    cum_probs <- exp(eta_matrix) / (1+exp(eta_matrix))
    # Initialize the matrix for category-specific probabilities
    p_vals <- matrix(nrow = n2_bar, ncol = num_categories)
    p_vals[, 1] <- cum_probs[, 1]
    for (j in 2:(num_categories - 1)) {
      p_vals[, j] <- cum_probs[, j] - cum_probs[, j - 1]
    }
    p_vals[, num_categories] <- 1 - cum_probs[, num_categories - 1]
  } else if (model_type == "Adjacent_Category") {


    # Calculate the cumulative sum for each row (from right to left)
    exp_cum_sum_matrix <- exp(t(apply(eta_matrix, 1, function(x) rev(cumsum(rev(x))))))


    p_vals <- matrix(nrow = n2_bar, ncol = num_categories)

    for (i in 1:n2_bar) {
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


    # Sum the Beta0 matrix and G matrix
    exp_sum_matrix <- exp(eta_matrix)

    p_vals <- matrix(nrow = n2_bar, ncol = num_categories)

    # Iterate over each row to compute probabilities
    for (i in 1:n2_bar) {
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

    # Sum the Beta0 matrix and G matrix
    exp_sum_matrix <- cbind(rep(1,n2_bar),exp(eta_matrix))
    # generate probabilities
    p_vals <- matrix(nrow = n2_bar, ncol = num_categories)
    for (j in 1:num_categories) {
      p_vals[, j] <- exp_sum_matrix[, j] / rowSums(exp_sum_matrix)
    }
    # p_vals[, 1] <- 1- rowSums(p_vals[, 2:(num_categories)])
  }















  # Compute the weights based on the computed probabilities
  dat0bis$wg0_1 <- sapply(1:nrow(dat0bis), function(i) {
    probs_i <- p_vals[i, ]
    x_vec <- integer(num_categories)
    x_vec[dat0bis[[resp]][i]] <- 1
    # Before the dmultinom line:
    if(any(!is.finite(probs_i)) || any(probs_i < 0) || all(probs_i == 0)) {
      print(paste("Row:", i))
      print(probs_i)
    }
    dmultinom(x = x_vec, size = 1, prob = probs_i)
  }) * dat0bis[["q"]]

  # Aggregate the weights
  dat0bis[, N1 := sum(wg0_1), by = id_]

  # Normalize the weights
  dat0bis$wg1 <- dat0bis$wg0_1 / dat0bis$N1

  # Order the data by the original order
  dat0bis <- dat0bis[order(dat0bis$ord),]

  # Return the computed weights
  return(dat0bis$wg1)
}
