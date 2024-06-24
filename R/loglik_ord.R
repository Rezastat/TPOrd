loglik_ord<-function(theta,q,formula,Y_,G_,Z_,dat,family,num_categories,model_type){
  datbis=merge(dat,q,by=c(G_,Z_))
  In<-t(fac2sparse(datbis$id_))#class.ind2.Mat(datbis$id_)

  n <- nrow(datbis)




  #########
  X0 = model.matrix(formula, datbis)
  extended_X0 <-  expand_model_matrix(X0, num_categories-1)





  if (model_type == "Stereotype_Regression"){
    phi_col <- rep(c(theta[1:(num_categories-2)],1),n)
    extended_X0[, num_categories:ncol(extended_X0)] <- extended_X0[, num_categories:ncol(extended_X0)] * phi_col
    eta_vals <- extended_X0 %*% theta[(num_categories-1):length(theta)]
  } else {
    eta_vals <- extended_X0 %*% theta
  }
  eta_matrix <- matrix(eta_vals,nrow= n,ncol = num_categories-1, byrow = T)

  if (model_type == "Proportional_Odds") {

    cum_probs <- exp(eta_matrix) / (1+exp(eta_matrix))
    # Initialize the matrix for category-specific probabilities
    p_vals <- matrix(nrow = n, ncol = num_categories)
    p_vals[, 1] <- cum_probs[, 1]
    for (j in 2:(num_categories - 1)) {
      p_vals[, j] <- cum_probs[, j] - cum_probs[, j - 1]
    }
    p_vals[, num_categories] <- 1 - cum_probs[, num_categories - 1]
  } else if (model_type == "Adjacent_Category") {


    # Calculate the cumulative sum for each row (from right to left)
    exp_cum_sum_matrix <- exp(t(apply(eta_matrix, 1, function(x) rev(cumsum(rev(x))))))


    p_vals <- matrix(nrow = n, ncol = num_categories)

    for (i in 1:(n)) {
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

    p_vals <- matrix(nrow = n, ncol = num_categories)

    # Iterate over each row to compute probabilities
    for (i in 1:(n)) {
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
    exp_sum_matrix <- cbind(rep(1,n),exp(eta_matrix))
    # generate probabilities
    p_vals <- matrix(nrow = n, ncol = num_categories)
    for (j in 1:num_categories) {
      p_vals[, j] <- exp_sum_matrix[, j] / rowSums(exp_sum_matrix)
    }
    # p_vals[, 1] <- 1- rowSums(p_vals[, 2:(num_categories)])
  }


  #########



  inner_arg <-   sapply(1:n, function(i) {
    probs_i <- p_vals[i, ]
    x_vec <- integer(num_categories)
    x_vec[datbis$Y[i]] <- 1
    dmultinom(x = x_vec, size = 1, prob = probs_i)
  }) *datbis$q*In

  ll <- sum( log( colSums( inner_arg ) ) )






  return(ll)
}
