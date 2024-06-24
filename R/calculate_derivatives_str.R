calculate_derivatives_str <- function(y_matrix, x_matrix, alpha, phi, beta, weight=NULL) {
  n <- nrow(y_matrix)
  h <- ncol(y_matrix)
  m <- ncol(x_matrix)

  # Create an empty matrix to store derivatives. Total columns = (h-2) + (h-1) + m
  derivative_matrix <- matrix(0, nrow = n, ncol = (h-2) + (h-1) + m)

  for (i in 1:n) {
    # Identify the category f for this observation
    f <- which(y_matrix[i,] == 1)

    # Common term used in all derivatives
    sum_beta_x <- sum(beta * x_matrix[i,])
    common_denominator <- sum(exp(alpha + phi * sum_beta_x))

    # Derivatives with respect to phi_f
    for (pf in 2:(h-1)) {
      if (pf == f) {
        derivative_matrix[i,  pf-1] <- sum_beta_x - exp(alpha[pf] + phi[pf] * sum_beta_x) * sum_beta_x / common_denominator
      } else{
        derivative_matrix[i,  pf-1] <-  - exp(alpha[pf] + phi[pf] * sum_beta_x) * sum_beta_x / common_denominator

      }
    }

    # Derivatives with respect to alpha_f
    for (af in 2:h) {
      if (af == f) {
        derivative_matrix[i, h-2 + af-1] <- 1 - exp(alpha[af] + phi[af] * sum_beta_x) / common_denominator
      } else {
        derivative_matrix[i, h-2 + af-1] <- - exp(alpha[af] + phi[af] * sum_beta_x) / common_denominator
      }
    }



    # Derivatives with respect to beta_g
    for (bg in 1:m) {
      derivative_matrix[i, (h-1) + (h-2) + bg] <- phi[f] * x_matrix[i, bg] - sum( exp(alpha + phi * sum_beta_x) * phi * x_matrix[i, bg]) / common_denominator
    }

  }
  if (length(weight)!=0){
    derivative_matrix <- weight * derivative_matrix
  }
  return(derivative_matrix)
}
