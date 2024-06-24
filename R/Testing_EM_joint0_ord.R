Testing_EM_joint0_ord <- function(theta,q,Betas_ids=NULL,G_,FIM.obj){
  if (is.null(Betas_ids)) {
    Betas_ids <- which(names(theta) %in% G_)}
  nbetas <- length(theta)
  phi <- c(theta, q[, "q"])
  n_phi <- length(phi)
  Q1 <- FIM.obj$FIM
  D <- rbind(diag(n_phi - 1), 0)
  F1 <- t(D) %*% Q1 %*% D
  Omega1 <- tryCatch(solve(F1), error = function(e) {
    ginv(F1)
  })
  IndU <- rep(F, n_phi)
  IndU[Betas_ids] <- T
  U <- FIM.obj$U
  U1 <- U[IndU]
  Ind <- rep(F, n_phi - 1)
  Ind[Betas_ids] <- T
  invF1_noInd <- tryCatch(solve(F1[!Ind, !Ind]), error = function(e) {
    ginv(F1[!Ind, !Ind])
  })
  V1 <- F1[Ind, Ind] - F1[Ind, !Ind] %*% invF1_noInd %*% F1[!Ind,
                                                            Ind]
  invOmega1_Ind <- tryCatch(solve(Omega1[Ind, Ind]), error = function(e) {
    ginv(Omega1[Ind, Ind])
  })
  W <- as.numeric(t(phi[IndU]) %*% invOmega1_Ind %*% phi[IndU])
  invV1 <- tryCatch(solve(V1), error = function(e) {
    ginv(V1)
  })
  S <- as.numeric(t(U1) %*% invV1 %*% U1)
  return(list(Var = diag(Omega1)[1:nbetas], W = W, S = S, df = length(Betas_ids)))
}
