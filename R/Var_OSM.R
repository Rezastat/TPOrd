Var_OSM <- function(obj_Ho, dat, formul){

  X1 <- model.matrix(formul, data = dat)
  x_matrix <- X1[,-1,drop=F]
  y1 <- dat[,"Y"]
  y_matrix <- as.matrix(model.matrix(~factor(y1)-1))


  alpha <- obj_Ho$alpha[, 1]
  phi <- obj_Ho$phi[, 1]
  beta <- obj_Ho$beta[, 1]

  I_in <- -calculate_hessian_str(y_matrix, x_matrix, alpha, phi, beta)


  return(diag(solve(I_in)))
}
