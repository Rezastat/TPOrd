optimTP_LM_ord  <-  function(formula,miscov,auxvar,strata,n,data,beta,p_gz,optimMeasure,K.idx=NULL,min.nk=NULL,logical.sub=NULL, num_categories,model_type){



  ### TODO: check if all unique values in data[,Z_] are in p_gz

  terms_dat <- unique( c(all.vars(formula), attr(terms(miscov), "term.labels") , attr(terms(auxvar), "term.labels"), attr(terms(strata), "term.labels") , attr(terms(auxvar), "term.labels")) )

  Y_ <- as.character(formula[[2]])
  G_ <- all.vars(miscov)
  Z_ <- all.vars(auxvar)
  S_ <- all.vars(strata)



  ## remove enties with zero or negative probabilities
  p_gz <- p_gz[p_gz$q>0,]

  N=NROW(data)
  rho=n/N
  if( rho<0.1 )warning("The phase 2 sample size (n2) is smaller than recommended. Please proceed with caution.")

  stratadf <- data.frame(xtabs(strata,data=data))
  #Xs <- model.matrix(strata, stratadf)
  stratadf[S_] <- lapply(stratadf[S_], function(x)as.numeric(as.character(x))) ## may need it for data.table
  g.vals <- unique(p_gz[,G_])
  z.vals <- unique(p_gz[,Z_])

  # if( !is.null(logical.sub) ){
  #   data <- data[logical.sub,]
  #   N <- NROW(data)
  # }

  dat_ <- cbind(id_=rep(seq(1,N,by=1), each=NROW(g.vals)),data[rep(seq_len(N), each=NROW(g.vals)), ],matrix(rep(t(g.vals),N), ncol=NCOL(g.vals), byrow=TRUE))
  names(dat_)[(ncol(dat_)-NCOL(g.vals)+1):ncol(dat_)] <- G_
  dat_ <- droplevels(dat_)
  dat_$ord_ <- 1:nrow(dat_)



  wgs <- omega1_g_ord(formula,model_type,dat_,beta,p_gz,G_,Z_) # P(G|Y,Z) and P(Y,Z)

  dat_$wg_ <- wgs$wg
  dat_ <- dat_[!is.na(dat_$wg_),]



  M1 <- obsIMf1_ord(formula,dat_,beta,p_gz,G_,Z_,by.id=FALSE,model_type=model_type, num_categories= num_categories)$I3
  M1 <- M1/N

  ####################

  results_list <- list()


  unique_values <- unique(dat_[, S_])


  for (i in seq_along(unique_values)) {


    df.by <- dat_[dat_[, S_] == unique_values[i], ]


    Nk <- length(unique(df.by$id_))
    IM <- obsIMf1_ord(formula, df.by, beta, p_gz, G_, Z_, model_type=model_type, num_categories=num_categories)


    result_df <- data.frame(df.by[1, S_, drop=FALSE], t(as.vector(IM$I3) / Nk))


    results_list[[i]] <- result_df
  }


  IM3_subs <- do.call(rbind, results_list)

  ########################



  results_list <- list()





  for (i in seq_along(unique_values)) {

    # Subset data for the current unique value
    df.by <- dat_[dat_[, S_] == unique_values[i], ]


    Nk <- length(unique(df.by$id_))
    IM <- obsIMf1_ord(formula, df.by, beta, p_gz, G_, Z_, model_type=model_type, num_categories=num_categories)


    result_df <- data.frame(df.by[1, S_, drop=FALSE], t(as.vector(IM$I2) / Nk))


    results_list[[i]] <- result_df
  }


  IM2_subs <- do.call(rbind, results_list)

  lenS_ <- length(S_)
  npars <- NCOL(M1)

  ### Indicator for the beta parameters of interest
  Ind <- rep(FALSE,npars-1)
  Betas_ids <- which(grepl(paste0(paste0("^",G_),collapse="|"),colnames(model.matrix(formula,dat_[1,]))))
  Ind[Betas_ids] <- TRUE

  if( optimMeasure=="Par-spec" & is.null(K.idx) ) stop("For a parameter-specific criterion K.idx must be provided.")
  ## asign these functions to the execution environment
  optMsuref=twoPhaseGAS:::.optMsure(optimMeasure);  environment(optMsuref)=environment()

  ### parameters for the constained optimization
  # if( !all.equal(pS_[S_],stratadf[,S_],attributes=FALSE,check.attributes = FALSE) ) stop("Strata groups in pS and stratadf don't match!")
  # ais <- (stratadf$Freq)#/sum(stratadf$Freq)
  # ais <- N*pS_$pS_#/sum(stratadf$Freq)
  npis <- nrow(stratadf)-1
  rho <- n
  upperLM <- as.numeric(stratadf$Freq) #rep(1,npis+1)
  # M1inv <- MASS::ginv(M1)
  # In <- diag(npars)
  D <- rbind(diag(npars-1),0)
  # Dinv <- MASS::ginv(D)
  # Dtinv  <- MASS::ginv((t(D)))
  ### result on how to calcuate the inverse of A+B
  # https://math.stackexchange.com/questions/17776/inverse-of-the-sum-of-matrices
  # (A+B)^{-1} = A^{-1}(I-(I+BA^{-1})^{-1}BA^{-1})


  ### The LM approach is going to get optimal values for nk/P(R=1,Z,K) as opposed to pr(R=1|Y,Z). In additio the  objective function is determined according to whether the design regression parameters are under the null or alternative

  if( all(beta[Ind]==0) ){
    ObjFun <- function(pis){ # pis = (upperLM*n/N)[-1]
      pis1 <- c(0,pis)
      pis1[1] <- (rho-sum(pis1))
      stratadf$prReq1YZ <- (pis1/N)/stratadf$Freq

      IM2b_subs <- as.data.frame(merge(as.data.table(IM2_subs), as.data.table(stratadf[,c(S_,"prReq1YZ")]), by=S_, sort = F))
      IM3b_subs <- as.data.frame(merge(as.data.table(IM3_subs), as.data.table(stratadf[,c(S_,"prReq1YZ")]), by=S_, sort = F))

      FIM <- M1 + matrix(colSums(IM2b_subs[,(lenS_+1):(npars^2+lenS_)]*IM2b_subs$prReq1YZ),ncol=npars) - matrix(colSums(IM3b_subs[,(lenS_+1):(npars^2+lenS_)]*IM3b_subs$prReq1YZ),ncol=npars)

      ## to deal with the constraint on the p_g's
      F1 <-  crossprod(D, FIM %*% D)

      ###  Score test variance
      invF1_noInd <- ginv(F1[!Ind,!Ind])
      V1 <- F1[Ind,Ind] - ( F1[Ind,!Ind] %*% invF1_noInd %*% F1[!Ind,Ind] )
      invV1 <- solve( V1 )

      return( optMsuref(invV1,K.idx) )
    }
  }else {
    ObjFun <- function(pis){ # pis = (upperLM*n/N)[-1]
      pis1 <- c(0,pis)
      pis1[1] <- (rho-sum(pis1))
      stratadf$prReq1YZ <- (pis1/N)/stratadf$Freq

      IM2b_subs <- as.data.frame(merge(as.data.table(IM2_subs), as.data.table(stratadf[,c(S_,"prReq1YZ")]), by=S_, sort = F))
      IM3b_subs <- as.data.frame(merge(as.data.table(IM3_subs), as.data.table(stratadf[,c(S_,"prReq1YZ")]), by=S_, sort = F))

      FIM <- M1 + matrix(colSums(IM2b_subs[,(lenS_+1):(npars^2+lenS_)]*IM2b_subs$prReq1YZ),ncol=npars) - matrix(colSums(IM3b_subs[,(lenS_+1):(npars^2+lenS_)]*IM3b_subs$prReq1YZ),ncol=npars)

      ## to deal with the constraint on the p_g's
      F1 <-  crossprod(D, FIM %*% D)

      # B <- matrix(colSums(IM2b_subs[,(lenS_+1):(npars^2+lenS_)]*IM2b_subs$prReq1YZ),ncol=npars) - matrix(colSums(IM3b_subs[,(lenS_+1):(npars^2+lenS_)]*IM3b_subs$prReq1YZ),ncol=npars)
      #
      #
      # FIMinv <- M1inv %*% (In - solve(In+B%*%M1inv)) %*% (B%*%M1inv)
      #
      # vcov = Dinv %*% FIMinv %*% Dtinv

      vcov <- tryCatch({
        solve(F1)
      }, error = function(e) {
        warning("Standard inversion failed due to singularity, using pseudoinverse instead.")
        ginv(F1)
      })

      return( optMsuref(vcov,K.idx) )

    }
  }


  if( is.null(min.nk) ){
    minLB.val <- 0 #ifelse(rho<0.05,rho*0.01,0.05)
    minLB <- rep(minLB.val,npis+1)
  }else{
    if( !length(min.nk) %in% c(1,npis+1) ) stop("The number of elements of min.nk is either 1 or ",npis+1)
    if( any(min.nk<0) )stop("The minimum stratum size must greater or equal than zero")
    if( any(min.nk>stratadf$Freq) )stop("The minimum stratum size must not be greater than the strata sample sizes determined by strata")
    minLB <- min.nk
  }
  initX0 <- (upperLM*n/N)[-1] #rep(rho,npis)
  initX1 <- (upperLM*n/N)  #rep(rho,npis+1)

  ObjFun_in <- function(pis){
    ui <- rbind(rep(-1,npis),rep(1,npis))
    ci <- c(-rho,0)
    return(as.numeric(ui%*%pis-ci))
  }

  ObjFun_eq <- function(pis1){
    return(sum(pis1)-rho)
  }

  ### Mind the warning: "For consistency with the rest of the package the inequality sign may be switched from >= to <= in a future nloptr version."

  ## To disable this warning and avoid printing do the following. However, new versions of the package may return an error when  this is enforced.
  options(nloptr.show.inequality.warning = FALSE)

  ### Note that the function can be changed when the change is implemented using
  # hin <- function(x) (-1)*f2(x, ...)  # NLOPT expects hin <= 0

  ### use this one for base, if an error then use the rest
  sol <- nloptr:::cobyla(initX0, fn = ObjFun, lower = minLB[-1], upper = upperLM[-1], hin = ObjFun_in, nl.info = FALSE, control = list(xtol_rel = 1e-6, maxeval = 10000))

  # return the option to the original state
  options(nloptr.show.inequality.warning = TRUE)

  if( !(sol$convergence>0 & sol$convergence<5) ){
    if( sol$convergence==5 ) message("The maximum number of evaluations in nloptr:::cobyla() has been reached. Trying dfoptim:::hjkb() now. You could also increase the number of evaluations via maxeval (although it is fairly large already.)")
    if( sol$convergence<0 ) message(paste0("Convergence in nloptr:::cobyla() was not achieved; error code ",sol$convergence,". Trying dfoptim:::hjkb() now."))

    ObjFunb <- function(pisA, k){ # pisA <- (upperLM*n/N)
      nstar <- sum(pisA)
      return( ObjFun(pisA[-1]) + k*(nstar-rho)^2 )
    }

    sol1 <- hjkb(initX1, fn=ObjFunb, lower=minLB, upper=upperLM, k=0.01, control=list(maxfeval=10000))
    sol1 <- tryCatch(hjkb(sol1$par, fn=ObjFunb, lower=minLB, upper=upperLM, k=100, control=list(maxfeval=10000)), error=function(e){ NULL })
    if( is.null(sol1) ){
      solB <- c(0,sol$par)
      solB[1] <- (rho-sum(solB))

      sol1 <- hjkb(solB, fn=ObjFunb, lower=minLB, upper=upperLM, k=100, control=list(maxfeval=10000))
    }

    if( sol1$convergence!=0 ){
      warning(paste0("Convergence in dfoptim:::hjkb() was not achieved. Error code ",sol1$convergence,". Returning last par value."))
    }

    stratadf$convergence <- c(sol1$convergence,rep(NA,nrow(stratadf)-1))
    pisopt <- sol1$par

  } else {
    stratadf$convergence <- c(sol$convergence,rep(NA,nrow(stratadf)-1))
    pisopt <- c(0,sol$par)
    pisopt[1] <- (rho-sum(pisopt))
  }

  stratadf$prR_cond_optim <- pisopt/sum(pisopt)

  return(stratadf)
}
