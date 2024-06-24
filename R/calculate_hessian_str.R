calculate_hessian_str <- function(y_matrix, x_matrix, alpha, phi, beta, weight=NULL) {
  n <- nrow(y_matrix)
  h <- ncol(y_matrix)
  m <- ncol(x_matrix)
  p <- (h-2) + (h-1) + m  # Total number of parameters

  # Initialize the total Hessian matrix
  total_hessian <- matrix(0, nrow = p, ncol = p)

  for (i in 1:n) {
    # Identify the category f for this observation
    f <- which(y_matrix[i,] == 1)

    # Calculate common terms
    sum_beta_x <- sum(beta * x_matrix[i,])
    exp_terms <- exp(alpha + phi * sum_beta_x)
    common_denominator <- sum(exp_terms)


    # Initialize the Hessian matrix for this observation
    hessian_matrix <- matrix(0, nrow = p, ncol = p)

    # Calculate second derivatives
    # Alpha derivatives
    for (af in 2:h) {

      # d2l/dalpha_f^2
      hessian_matrix[(h-2+af-1), (h-2+af-1)] <-  -(exp(alpha[af]+phi[af]*sum_beta_x)*common_denominator-exp(alpha[af]+phi[af]*sum_beta_x)^2)/common_denominator^2


    }
    # Alpha-Alpha
    for (af in 2:h){

      if((af+1) <= h){
        for(ad in (af+1):h){
          # d2l/dalpha_f dalpha_d
          hessian_matrix[(h-2+af-1), (h-2+ad-1)] <- (exp(alpha[af]+phi[af]*sum_beta_x) * exp(alpha[ad]+phi[ad]*sum_beta_x)) / common_denominator^2
        }
      }

    }
    # Alpha-Beta
    for (af in 2:h) {

      for (bg in 1:m) {
        # d2l/dalpha_f dbeta_g
        hessian_matrix[(h-2+af-1), (h-2+h-1 + bg)] <- -(exp(alpha[af]+phi[af]*sum_beta_x) * phi[af] * x_matrix[i, bg] * common_denominator - (sum(exp_terms * phi * x_matrix[i, bg]))*exp(alpha[af]+phi[af]*sum_beta_x) )/ common_denominator^2

      }

    }
    # Alpha-Phi
    for (af in 2:h) {


      for (pd in 2:(h-1)){
        if (af == pd) {
          # d2l/dalpha_f dphi_f
          hessian_matrix[(pd-1), (h-2 + af-1)] <- -(exp(alpha[af]+phi[af]*sum_beta_x) * sum_beta_x * common_denominator - exp(alpha[af]+phi[af]*sum_beta_x)^2 * sum_beta_x) / common_denominator^2
        }else {
          hessian_matrix[(pd-1), (h-2 + af-1)] <- (exp(alpha[af]+phi[af]*sum_beta_x)* exp(alpha[pd]+phi[pd]*sum_beta_x)* sum_beta_x)/common_denominator^2
        }
      }

    }



    # Phi^2
    for (pf in 2:(h-1)) {
      # d2l/dphi_f^2

      hessian_matrix[(pf-1), (pf-1)] <- -(sum_beta_x^2 * exp(alpha[pf]+phi[pf]*sum_beta_x) * common_denominator - (sum_beta_x * exp(alpha[pf]+phi[pf]*sum_beta_x))^2) / common_denominator^2


    }
    # Phi-Beta
    for (pf in 2:(h-1)) {

      for (bg in 1:m){
        if (pf==f){
          hessian_matrix[(pf-1),(h-2+h-1+bg)] <- x_matrix[i,bg] - ((exp(alpha[pf]+phi[pf]*sum_beta_x)*sum_beta_x *phi[pf]*x_matrix[i,bg]+exp(alpha[pf]+phi[pf]*sum_beta_x)*x_matrix[i,bg])*common_denominator-exp(alpha[pf]+phi[pf]*sum_beta_x)*sum_beta_x*sum(exp_terms*phi*x_matrix[i,bg]))/common_denominator^2
        } else {
          hessian_matrix[(pf-1),(h-2+h-1+bg)] <-  - ((exp(alpha[pf]+phi[pf]*sum_beta_x)*sum_beta_x *phi[pf]*x_matrix[i,bg]+exp(alpha[pf]+phi[pf]*sum_beta_x)*x_matrix[i,bg])*common_denominator-exp(alpha[pf]+phi[pf]*sum_beta_x)*sum_beta_x*sum(exp_terms*phi*x_matrix[i,bg]))/common_denominator^2
        }
      }

    }

    # Phi-Phi
    for (pf in 2:(h-1)){

      if (pf +1 <= h-1){
        for (pd in (pf+1):(h-1)){
          hessian_matrix[(pf-1),(pd-1)] <- (exp(alpha[pf]+phi[pf]*sum_beta_x)*sum_beta_x^2* exp(alpha[pd]+phi[pd]* sum_beta_x))/common_denominator^2
        }

      }
    }
    # Beta^2
    for (bg in 1:m) {

      hessian_matrix[(h-2+h-1+ bg), (h-2+h-1+ bg)] <- -(sum(exp_terms * phi^2 * x_matrix[i, bg]^2)*common_denominator-sum(exp_terms * phi * x_matrix[i, bg])^2) / common_denominator^2


    }

    # Beta-Beta
    for (bg in 1:m) {
      if (bg + 1 <= m) {
        for (bc in (bg +1):m){
          hessian_matrix[(h-2+h-1+ bg), (h-2+h-1+ bc)] <- - (sum(phi^2* x_matrix[i,bg]*x_matrix[i,bc]*exp_terms) * common_denominator-sum(phi* x_matrix[i,bg]*exp_terms)*sum(phi* x_matrix[i,bc]*exp_terms))/common_denominator^2
        }
      }
    }







    # Add the Hessian of this observation to the total Hessian
    if(length(weight)==0){
      total_hessian <- total_hessian + hessian_matrix
    } else {
      total_hessian <- total_hessian + weight[i] * hessian_matrix
    }

  }

  return(total_hessian+t(total_hessian)-diag(diag(total_hessian)))
}
