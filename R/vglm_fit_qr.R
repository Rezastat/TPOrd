vglm_fit_qr <- function (x, y, w = rep_len(1, nrow(x)), X.vlm.arg = NULL, Xm2 = NULL,
                         Ym2 = NULL, etastart = NULL, mustart = NULL, coefstart = NULL,
                         offset = 0, family, control = vglm.control(), qr.arg = FALSE,
                         constraints = NULL, extra = NULL, Terms = Terms, function.name = "vglm",
                         ...){
  suppressMessages(attach(asNamespace("VGAM")))

  if (is.null(criterion <- control$criterion)){
    criterion <- "coefficients"}
  eff.n <- nrow(x)
  specialCM <- NULL
  post <- list()
  check.rank <- control$Check.rank
  nonparametric <- FALSE
  epsilon <- control$epsilon
  maxit <- control$maxit
  save.weights <- control$save.weights
  trace <- control$trace
  orig.stepsize <- control$stepsize
  minimize.criterion <- control$min.criterion
  history <- NULL
  fv <- NULL
  n <- nrow(x)
  stepsize <- orig.stepsize
  old.coeffs <- coefstart
  intercept.only <- ncol(x) == 1 && colnames(x) == "(Intercept)"
  y.names <- predictors.names <- NULL
  n.save <- n
  if (length(slot(family, "initialize"))) {
    eval(slot(family, "initialize"))}
  if (length(etastart)) {
    eta <- etastart
    mu <- if (length(mustart)) {
      mustart}
    else{ slot(family, "linkinv")(eta, extra = extra)}
  }
  if (length(mustart)) {
    mu <- mustart
    if (length(body(slot(family, "linkfun")))) {
      eta <- slot(family, "linkfun")(mu, extra = extra)
    }
    else {
      warning("argument 'mustart' assigned a value ", "but there is no 'linkfun' slot to use it")
    }
  }
  validparams <- validfitted <- TRUE
  if (length(body(slot(family, "validparams")))) {
    validparams <- slot(family, "validparams")(eta, y = y,
                                               extra = extra)}
  if (length(body(slot(family, "validfitted")))) {
    validfitted <- slot(family, "validfitted")(mu, y = y,
                                               extra = extra)}
  if (!(validparams && validfitted)) {
    stop("could not obtain valid initial values. ", "Try using 'etastart', 'coefstart' or 'mustart', else ",
         "family-specific arguments such as 'imethod'.")}
  M <- NCOL(eta)
  if (length(slot(family, "constraints"))){
    eval(slot(family, "constraints"))}
  Hlist <- process.constraints(constraints, x = x, M = M, specialCM = specialCM,
                               Check.cm.rank = control$Check.cm.rank)
  ncolHlist <- unlist(lapply(Hlist, ncol))
  X.vlm.save <- if (length(X.vlm.arg)) {
    X.vlm.arg
  }else {
    X.vlm.save <-  lm2vlm.model.matrix(x, Hlist, xij = control$xij, Xm2 = Xm2)
    if (length(coefstart)) { # this is 0 so doesn't run
      eta <- if (ncol(X.vlm.save) > 1) {
        matrix(X.vlm.save %*% coefstart, n, M, byrow = TRUE) +
          offset
      }  else {
        matrix(X.vlm.save * coefstart, n, M, byrow = TRUE) +
          offset
      }
      if (M == 1) {
        eta <- c(eta)}
      mu <- slot(family, "linkinv")(eta, extra = extra)
    }
    if (criterion != "coefficients") {
      tfun <- slot(family, criterion)
    } # !!! deviance
    iter <- 1
    new.crit <- switch(criterion, coefficients = 1, tfun(mu = mu,
                                                         y = y, w = w, res = FALSE, eta = eta, extra = extra))
    # important stuff starts here
    deriv.mu <- eval(slot(family, "deriv"))
    wz <- eval(slot(family, "weight"))
    # return(wz) #####
    if (control$checkwz) {
      wz <- checkwz(wz, M, trace = trace, wzepsilon = control$wzepsilon)}
    U <- vchol(wz, M = M, n = n, silent = !trace)
    tvfor <- vforsub(U, as.matrix(deriv.mu), M = M, n = n)
    z <- eta + vbacksub(U, tvfor, M = M, n = n) - offset
    one.more <- TRUE
    nrow.X.vlm <- nrow(X.vlm.save)
    ncol.X.vlm <- ncol(X.vlm.save)
    if (nrow.X.vlm < ncol.X.vlm) {
      stop("There are ", ncol.X.vlm, " parameters but only ",
           nrow.X.vlm, " observations")}
    #while (one.more) {
    ##!!!!qr
    tfit <- VGAM:::vlm.wfit(xmat = X.vlm.save, zmat = z, Hlist = NULL,
                            U = U, matrix.out = FALSE, is.vlmX = TRUE, qr = qr.arg,
                            xij = NULL)}
  # detach("package:VGAM")
  return(list(tfit=tfit,deriv.mu=deriv.mu))}
