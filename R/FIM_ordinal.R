FIM_ordinal <- function(theta,qG,formula,dat,fit=NULL,model_type, n_second = n_second, num_categories){
  X1 <- model.matrix(formula,dat)
  y1 <- dat[,"Y"]
  y_matrix <- as.matrix(model.matrix(~factor(y1)-1))
  nq<-nrow(qG); nbetas<-length(theta)
  phi<-c(theta,log(qG[,"q"])) #send the log(q_j) here to deal with the numerical issue on the hessian computation
  n_phi<-length(phi)

  Qd1=sapply(dat$k_,function(k){1/qG[qG[,"k_"]==k,"q"]})
  Qd2=sapply(dat$k_,function(k){-1/qG[qG[,"k_"]==k,"q"]^2})
  Qmat = as.matrix(t(as.matrix(Matrix::fac2sparse(as.factor(dat$k_)))))#class.ind2(dat$k_)
  Qmat.d1=Qmat*Qd1
  Psi=dat$wg_

  #######
  if (model_type == "Stereotype_Regression"){
    y1 <- dat[,"Y"]
    y_matrix <- as.matrix(model.matrix(~factor(y1)-1))

    x_matrix <- X1[,-1,drop=F]
    alpha <- c(0,(theta[(num_categories-1):(num_categories-1+num_categories-2)]))
    phi <- c(0,(theta[1:(num_categories-2)]),1)
    beta <- theta[(num_categories-1+num_categories-1):length(theta)]

    # Compute the derivatives
    I1_0 <- -calculate_hessian_str(y_matrix, x_matrix, alpha, phi, beta, weight = Psi)
  } else {

    ## second derivatives (i.e. hessian of log P(Y_i|g_j) + log q_j or jacobian of the gradient)

    if(!inherits(fit[["tfit"]][["qr"]], "qr")){

      I1_0<-crossprod(fit[["tfit"]][["qr"]])

    } else{
      I1_0<-crossprod(qr.R(fit[["tfit"]][["qr"]]))#/disp # upper left matrix with the fisher information matrix for betas

    }
  }


  #######





  I2_0 <-diag(colSums(Psi * Qmat * Qd2)) #second derivatives for an i,jth observation

  I1<-rbind(cbind(I1_0, matrix(0, nrow=nrow(I1_0), ncol=ncol(I2_0))),
            cbind(matrix(0, nrow=nrow(I2_0), ncol=ncol(I1_0)), -I2_0))

  ###### second matrix ###########
  if (model_type == "Stereotype_Regression"){
    sqrt_wl <- calculate_derivatives_str(y_matrix, x_matrix, alpha, phi, beta, weight = sqrt(Psi))
    pmat <- Qmat.d1*sqrt(Psi)
    sqrt_wl<- cbind(sqrt_wl,pmat)
    sum_matrix2 <- matrix(0, nrow = ncol(sqrt_wl), ncol = ncol(sqrt_wl))

    for (i in 1:nrow(sqrt_wl)) {

      outer_product <- sqrt_wl[i,] %*% t(sqrt_wl[i,])

      sum_matrix2 <- sum_matrix2 + outer_product
    }


  } else {
    sqrt_wl_extended<- dloglikelihood_dbeta(x_design=X1,y_matrix=y_matrix,thetas=theta,num_categories=num_categories,model_type,w=sqrt(Psi))
    sqrt_wl <- shrink_matrix(sqrt_wl_extended,num_categories)
    pmat <- Qmat.d1*sqrt(Psi)
    sqrt_wl<- cbind(sqrt_wl,pmat)
    sum_matrix2 <- matrix(0, nrow = ncol(sqrt_wl), ncol = ncol(sqrt_wl))


    for (i in 1:nrow(sqrt_wl)) {

      outer_product <- sqrt_wl[i,] %*% t(sqrt_wl[i,])

      sum_matrix2 <- sum_matrix2 + outer_product
    }

  }


  ########## third matrix #######
  if (model_type == "Stereotype_Regression"){
    wl <- calculate_derivatives_str(y_matrix, x_matrix, alpha, phi, beta, weight = Psi)
  } else {
    wl_extended<- dloglikelihood_dbeta(x_design=X1,y_matrix=y_matrix,thetas=theta,num_categories=num_categories,model_type,w=Psi)
    wl <- shrink_matrix(wl_extended,num_categories)
  }


  pmat <- Qmat.d1*Psi
  wl <- cbind(wl, pmat)
  U <- colSums(wl)
  wl_summed <- sum_kth_rows(wl,(n_second+1), max(unique(dat$G))) ###### change 2001 !!!


  sum_matrix3 <- matrix(0, nrow = ncol(wl_summed), ncol = ncol(wl_summed))
  for (i in 1:nrow(wl_summed)) {

    outer_product3 <- wl_summed[i,] %*% t(wl_summed[i,])

    sum_matrix3 <- sum_matrix3 + outer_product3
  }

  ###############################


  obsFIM=I1-sum_matrix2+sum_matrix3
  return(list(FIM=obsFIM,I1=I1,I2=sum_matrix2,I3=sum_matrix3,U=U))
}
