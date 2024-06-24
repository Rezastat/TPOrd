fitnessTP_ord  <-  function(obsIM.R,Rj,optimMeasure,K.idx=NULL)  {

  FIM <- obsIM.R$IM1 + matrix(colSums(obsIM.R$IM2_id*Rj),ncol=obsIM.R$npars) - matrix(colSums(obsIM.R$IM3_id*Rj),ncol=obsIM.R$npars)

  D <- rbind(diag(obsIM.R$npars-1),0)
  F1 <-  crossprod(D, FIM %*% D)

  ### Function that determines the optimality measure (D-, A- or parameter-specfic for now)
  ### if optimMeasure=="Par-spec" then K.idx MUST be not null
  if( optimMeasure=="Par-spec" & is.null(K.idx) ) stop("For a parameter-specific criterion K.idx must be provided.")

  ## asign these functions to the execution environment
  optMsuref=twoPhaseGAS:::.optMsure(optimMeasure)

  if( obsIM.R$Null ){
    Ind <- obsIM.R$Ind

    invF1_noInd <- ginv(F1[!Ind,!Ind])
    V1 <- F1[Ind,Ind] - F1[Ind,!Ind] %*% invF1_noInd %*% F1[!Ind,Ind]
    invV1 <- solve( V1 )

    return( optMsuref(invV1,K.idx) )

  }else {
    vcov <- tryCatch({
      solve(F1)
    }, error = function(e) {
      warning("Standard inversion failed; using pseudoinverse")
      ginv(F1)
    })

    return(optMsuref(vcov, K.idx))
  }
}

