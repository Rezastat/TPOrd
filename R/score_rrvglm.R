score_rrvglm <- function(rrvglm_obj, rrvglm_obj_Ho, dat, num_categories){

  X1 <- rrvglm_obj@x
  coefs <- summary(rrvglm_obj_Ho)@coef3[,1]
  coefs <- c(coefs[(num_categories-2):1],coefs[(num_categories-1+num_categories-2):(num_categories-1)],0, coefs[(num_categories-1+num_categories-1):length(coefs)])

  y1 <- dat[,"Y"]
  y_matrix <- as.matrix(model.matrix(~factor(y1)-1))

  x_matrix <- X1[,-1,drop=F]
  alpha <- c(0,(coefs[(num_categories-1):(num_categories-1+num_categories-2)]))
  phi <- c(0,(coefs[1:(num_categories-2)]),1)
  beta <- coefs[(num_categories-1+num_categories-1):length(coefs)]

  I_in <- -calculate_hessian_str(y_matrix, x_matrix, alpha, phi, beta)

  S <- colSums(calculate_derivatives_str(y_matrix, x_matrix, alpha, phi, beta))
  STest <- t(S) %*% solve(I_in) %*% S
  return(STest)
}
