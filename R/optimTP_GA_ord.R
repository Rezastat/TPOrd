optimTP_GA_ord <- function(ncores, formula, miscov, auxvar, model_type, n, data, beta, p_gz, ga.popsize, ga.propelit, ga.proptourney, ga.ngen, ga.mutrate, ga.initpop = NULL, optimMeasure, K.idx = NULL, seed = 1, verbose = 0) {

  ptm0 <- Sys.time()
  terms_dat <- unique(c(all.vars(formula), attr(terms(miscov), "term.labels"), attr(terms(auxvar), "term.labels")))

  Y_ <- as.character(formula[[2]])
  G_ <- all.vars(miscov)
  Z_ <- all.vars(auxvar)

  ## Remove entries with zero or negative probabilities
  p_gz <- p_gz[p_gz$q > 0,]

  ### Check if names in G_ and Z_ are in the joint distribution data.frame p_gz
  g.vals <- unique(p_gz[, G_])
  N <- NROW(data)

  ### Function that determines the optimality measure (D-, A- or parameter-specific for now)
  ### If optimMeasure == "Par-spec" then K.idx MUST be not null
  if (optimMeasure == "Par-spec" & is.null(K.idx)) stop("For a parameter-specific criterion K.idx must be provided.")

  ## Assign these functions to the execution environment
  omega1_g_ord = omega1_g_ord
  optMsuref = twoPhaseGAS:::.optMsure(optimMeasure)

  envup = environment(omega1_g_ord) = environment(optMsuref) = environment()

  uniqZ <- unique(p_gz[, Z_])
####### problem here:
  uniqZ <- uniqZ[order(uniqZ)]
  # uniqZ <- uniqZ[order(uniqZ$Z), ]


  ### Calculate FIMs for all samples (hopefully, a faster way than before)
  dat_ <- cbind(id_ = rep(seq(1, N, by = 1), each = NROW(g.vals)), data[rep(seq_len(N), each = NROW(g.vals)), ], matrix(rep(t(g.vals), N), ncol = NCOL(g.vals), byrow = TRUE))
  names(dat_)[(ncol(dat_) - NCOL(g.vals) + 1):ncol(dat_)] <- G_
  dat_ <- droplevels(dat_)

  wgs <- omega1_g_ord(formula, model_type, dat_, beta, p_gz, G_, Z_)

  dat_$wg_ <- wgs$wg
  dat_ <- dat_[!is.na(dat_$wg_),]

  IM1 <- obsIMf1_ord(formula, dat_, beta, p_gz, G_, Z_, model_type = model_type, num_categories = num_categories)$I3 / N

  obsIM23_id <- obsIMf1_ord(formula, dat_, beta, p_gz, G_, Z_, model_type = model_type, num_categories = num_categories, by.id = TRUE)

  npars <- NCOL(IM1)

  IM2_id <- obsIM23_id$I2 / N
  IM3_id <- obsIM23_id$I3 / N

  Ind <- rep(FALSE, npars - 1)
  Betas_ids <- which(grepl(paste0(paste0("^", G_), collapse = "|"), colnames(model.matrix(formula, dat_[1, ]))))
  Ind[Betas_ids] <- TRUE
  D <- rbind(diag(npars - 1), 0)

  if (all(beta[Ind] == 0)) {
    Assess_fitness <- function(R) {
      FIM <- IM1 + matrix(colSums(IM2_id * R), ncol = npars) - matrix(colSums(IM3_id * R), ncol = npars)
      F1 <- t(D) %*% FIM %*% D

      invF1_noInd <- ginv(F1[!Ind, !Ind])
      V1 <- F1[Ind, Ind] - F1[Ind, !Ind] %*% invF1_noInd %*% F1[!Ind, Ind]
      invV1 <- tryCatch({
        solve(V1)
      }, error = function(e) {
        warning("Standard inversion failed; using pseudoinverse")
        ginv(V1)
      })

      return(optMsuref(invV1, K.idx))
    }
  } else {
    Assess_fitness <- function(R) {
      FIM <- IM1 + matrix(colSums(IM2_id * R), ncol = npars) - matrix(colSums(IM3_id * R), ncol = npars)
      F1 <- t(D) %*% FIM %*% D

      vcov <- tryCatch({
        solve(F1)
      }, error = function(e) {
        warning("Standard inversion failed; using pseudoinverse")
        ginv(F1)
      })
      return(optMsuref(vcov, K.idx))
    }
  }

  R0s <- rep(0, N)
  ObjFun <- function(v) {
    R <- R0s
    R[v] <- 1
    Assess_fitness(R)
  }

  time.elapsed0 <- as.numeric(difftime(Sys.time(), ptm0, units = "secs"))

  ptm <- Sys.time()

  out <- twoPhaseGAS:::kofnGA.mc(verbose = verbose, ncores, n = N, k = n, OF = ObjFun, popsize = ga.popsize,
                                 keepbest = ceiling(ga.popsize * ga.propelit),
                                 tourneysize = max(ceiling(ga.popsize * ga.proptourney), 2),
                                 ngen = ga.ngen, mutprob = ga.mutrate, initpop = ga.initpop,
                                 varsExport = c("Assess_fitness", "optMsuref", "K.idx", "IM1", "IM2_id", "IM3_id", "D", "R0s"), envup = envup)

  time.elapsed <- as.numeric(difftime(Sys.time(), ptm, units = "secs"))
  out$time.elapsed <- c(preproc = time.elapsed0, GA = time.elapsed)

  return(out)
}

