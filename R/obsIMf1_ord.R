obsIMf1_ord <- function(formula,dat,beta,p_gz,G_,Z_,by.id=FALSE,model_type,num_categories){
  Y <- dat[,as.character(formula[[2]])]
  X <- model.matrix(formula, dat)



  Omega <- dat$wg_


  # qG <- aggregate(as.formula(paste0("q~",paste(G_,collapse="+"))),p_gz,FUN=sum); names(qG)[which(names(qG)=="q")]<-"qG"
  qG <- p_gz; names(qG)[which(names(qG)=="q")]<-"qG"

  qG$k_ <- 1:nrow(qG)

  # dat$ord_ <- 1:nrow(dat)
  # # dat <- as.data.frame(merge(as.data.table(dat),as.data.table(qG[,c(G_,"qG","k_")]),by=G_, sort = F)) # merge(dat,qG[,c(G_,"qG","k_")],by=G_) ## this one may be sped up
  # dat <- as.data.frame(merge(as.data.table(dat),as.data.table(qG[,c(G_,Z_,"qG","k_")]),by=c(G_,Z_), sort = F)) # merge(dat,qG[,c(G_,"qG","k_")],by=G_) ## this one may be sped up
  # dat <- dat[order(dat$ord_),]
  # ## Need the 3 lines below to always have the same number of columns in Qmat irrespective of the subset
  # Qd1 <- 1/dat$qG
  # Qmat0 <- crossprod(sapply(dat$k_,function(x)1*x==qG$k_), diag(nrow(qG)))
  # Qmat.d1bis <- Qmat0*Qd1

  Qmat <- tcrossprod(matrix(rep(1,nrow(dat)),ncol=1),1/qG$qG)
  matchrows <- if( all(c(G_,Z_)%in%names(dat)) & all(c(G_,Z_)%in%names(qG)) ){
    match(interaction(dat[,c(G_,Z_)]), interaction(qG[,c(G_,Z_)]))
  }else if( Z_ %in% names(qG) & Z_ %in% names(dat) ){
    match(dat[,c(Z_)],qG[,c(Z_)])
  }else if( G_ %in% names(qG) & G_ %in% names(dat)){
    match(dat[,c(G_)],qG[,c(G_)])
  }else stop("No matching columns b/t p_gz and dat")

  Qmat.d1 <- Qmat*t(sapply(matchrows,function(s){1*s==1:nrow(qG)}))
  # Qmat.d1b <- crossprod(sapply(dat[,Z_],function(z){1*(p_gz[,Z_]==z)/p_gz[,"q"]}), diag(nrow(p_gz)))
  colnames(Qmat.d1) <- 1:nrow(qG)

  ## Original way assuming all the values of qG are in dat$k_
  # Qmat <- as.matrix(t(fac2sparse(dat$k_)))   #class.ind2(dat$k_)
  # Qmat.d1 <- Qmat*Qd1

  X <- X[, , drop = FALSE]

  if (model_type == "Stereotype_Regression"){
    y1 <- dat[,"Y"]


    if (length(unique(dat[, "Y"]))==1){
      y_matrix <- matrix(0, nrow = nrow(dat), ncol = num_categories)
      column_index <- unique(dat[, "Y"])
      y_matrix[, column_index] <- 1

    } else{
      y_matrix <- as.matrix(model.matrix(~factor(y1)-1))

    }


    x_matrix <- X[,-1,drop=F]
    alpha <- c(0,(beta[(num_categories-1):(num_categories-1+num_categories-2)]))
    phi <- c(0,(beta[1:(num_categories-2)]),1)
    betas <- beta[(num_categories-1+num_categories-1):length(beta)]
  } else{
    y1 <- dat[,"Y"]
    if (length(unique(dat[, "Y"]))==1){
      y_matrix <- matrix(0, nrow = nrow(dat), ncol = num_categories)
      column_index <- unique(dat[, "Y"])
      y_matrix[, column_index] <- 1

    } else{
      y_matrix <- as.matrix(model.matrix(~factor(y1)-1))

    }
  }

  if (model_type == "Stereotype_Regression"){
    first_Deriv <- calculate_derivatives_str(y_matrix, x_matrix, alpha, phi, betas)
    U <- cbind(first_Deriv,Qmat.d1)

    # for (i in 1:nrow(sqrt_wl)) {
    #
    #   outer_product <- sqrt_wl[i,] %*% t(sqrt_wl[i,])
    #
    #   sum_matrix2 <- sum_matrix2 + outer_product
    # }


  } else {

    first_Deriv<- dloglikelihood_dbeta(x_design=X,y_matrix=y_matrix,thetas=beta,num_categories=num_categories,model_type)
    first_Deriv <- sum_kth_rows(first_Deriv,1,(num_categories-1))
    U <- cbind(first_Deriv,Qmat.d1)

  }




  Ind_id <- fac2sparse(dat$id_)
  EU <- Ind_id %*% (Omega*U)

  if( !by.id ){

    obsI2 <- crossprod( sqrt(Omega)*U,  sqrt(Omega)*U )


    obsI3 <- as.matrix(crossprod( EU,  EU ))


  } else{
    rwouter <- function(mat){t(apply( mat, 1, function(x) x%o%x))}
    ## calculate second matrix, i.e. expected value of cross score functions
    obsI2 <- as.matrix( Ind_id %*% rwouter( U*sqrt(Omega) ) )
    ## Third matrix, i.e. outer product of expected score
    obsI3 <- rwouter( EU )
  }

  return(list(I2=obsI2,I3=obsI3))
}
