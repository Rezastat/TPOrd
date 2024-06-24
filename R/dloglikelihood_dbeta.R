dloglikelihood_dbeta <- function(x_design,y_matrix,thetas,num_categories,model_type,w=NULL){

  nrow_x <- nrow(x_design)
  ex_design_matrix <- expand_model_matrix(x_design,num_cat = num_categories-1)
  eta <- ex_design_matrix %*% thetas
  eta_matrix <- matrix(eta,nrow=nrow_x,ncol = num_categories-1, byrow = T)
  M <- NCOL(eta_matrix)
  # y1 <- dat[,"Y"]
  # y_matrix <- as.matrix(model.matrix(~factor(y1)-1))


  if (model_type == "Proportional_Odds"){
    cump <- cbind(eta2theta(eta_matrix, "logitlink", earg = list(
      theta =c() , bvalue = NULL, inverse = FALSE, deriv = 0,
      short = TRUE, tag = FALSE)), 1)
    mu <- cbind(cump[, 1], tapplymat1(cump, "diff"))

    mu.use <- pmax(mu, .Machine$double.eps * 1)

    dcump.deta <- dtheta.deta(cump, "logitlink", earg = list(
      theta = c(), bvalue = NULL, inverse = FALSE, deriv = 0,
      short = TRUE, tag = FALSE))[,-num_categories]



    if (length(w)==0){
      w<- rep(1, nrow_x)
    }




    deriv.answer <- w * dcump.deta * (y_matrix[, -(M +1), drop = FALSE]/mu.use[, -(M + 1), drop = FALSE] - y_matrix[, -1, drop = FALSE]/mu.use[, -1, drop = FALSE])
  } else if (model_type == "Adjacent_Category"){
    zeta <- eta2theta(eta_matrix, "loglink", list(theta = c(), bvalue = NULL,
                                                  inverse = FALSE, deriv = 0, short = TRUE, tag = FALSE))
    dzeta.deta <- dtheta.deta(zeta, "loglink", earg = list(theta = c(),
                                                           bvalue = NULL, inverse = FALSE, deriv = 0, short = TRUE,
                                                           tag = FALSE))
    d1 <- VGAM:::acat.deriv(zeta, M = M, n = nrow_x, reverse = TRUE)
    score <- attr(d1, "gradient")/d1
    if (length(w)==0){
      w<- rep(1, nrow_x)
    }
    deriv.answer <- if (TRUE) {
      cumy <- tapplymat1(y_matrix, "cumsum")
      c(w) * dzeta.deta * (cumy[, 1:M]/zeta - score)
    }

  } else if (model_type == "Stopping_Ratio") {


    mymat <- tapplymat1(y_matrix[, ncol(y_matrix):1, drop = FALSE], "cumsum")[, ncol(y_matrix):1, drop = FALSE]


    dj <- eta2theta(eta_matrix, "logitlink", earg = list(theta = c(),
                                                         bvalue = NULL, inverse = FALSE, deriv = 0, short = TRUE,
                                                         tag = FALSE))
    if (length(w)==0){
      w<- rep(1, nrow_x)
    }
    deriv.answer <- c(w) * (y_matrix[, -ncol(y_matrix), drop = FALSE]/dj - mymat[,
                                                                                 -1, drop = FALSE]/(1 - dj)) * dtheta.deta(dj, "logitlink",
                                                                                                                           earg = list(theta = c() , bvalue = NULL, inverse = FALSE,
                                                                                                                                       deriv = 0, short = TRUE, tag = FALSE))



  }

  expanded_dtheta_deta2<- as.vector(t(deriv.answer))
  dll_dbeta <- ex_design_matrix * expanded_dtheta_deta2
  return(dll_dbeta)
}
