vglm_qr <- function (formula, family = stop("argument 'family' needs to be assigned"),
                     data = list(), weights = NULL, subset = NULL, na.action = na.fail,
                     etastart = NULL, mustart = NULL, coefstart = NULL, control = vglm.control(...),
                     offset = NULL, method = "vglm.fit", model = FALSE, x.arg = TRUE,
                     y.arg = TRUE, contrasts = NULL, constraints = NULL, extra = list(),
                     form2 = NULL, qr.arg = TRUE, smart = TRUE, ...)
{
  dataname <- as.character(substitute(data))
  function.name <- "vglm"
  ocall <- match.call()
  if (smart)
    setup.smart("write")
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  parent.frame()
  switch(method, model.frame = return(mf), vglm.fit = 1, stop("invalid 'method': ",
                                                              method))
  mt <- attr(mf, "terms")
  xlev <- .getXlevels(mt, mf)
  y <- model.response(mf, "any")
  x <- if (!is.empty.model(mt)) {
    model.matrix(mt, mf, contrasts)
  }else {matrix(, NROW(y), 0)}
  oasgn <- attr(x, "assign")
  attr(x, "assign") <- VGAM:::attrassigndefault(x, mt)
  attr(x, "orig.assign.lm") <- oasgn
  if (!is.null(form2)) {
    if (!is.null(subset))
      stop("argument 'subset' cannot be used when ", "argument 'form2' is used")
    retlist <- shadowvglm(formula = form2, family = family,
                          data = data, na.action = na.action, control = vglm.control(...),
                          method = method, model = model, x.arg = x.arg, y.arg = y.arg,
                          contrasts = contrasts, constraints = constraints,
                          extra = extra, qr.arg = qr.arg)
    Ym2 <- retlist$Ym2
    Xm2 <- retlist$Xm2
    if (length(Ym2)) {
      if (NROW(Ym2) != NROW(y))
        stop("number of rows of 'y' and 'Ym2' are unequal")
    }
    if (length(Xm2)) {
      if (NROW(Xm2) != NROW(x))
        stop("number of rows of 'x' and 'Xm2' are unequal")
    }
  }
  else {
    Xm2 <- Ym2 <- NULL
  }
  offset <- model.offset(mf)
  if (is.null(offset))
    offset <- 0
  w <- model.weights(mf)
  if (!length(w)) {
    w <- rep_len(1, nrow(mf))
  }
  else if (NCOL(w) == 1 && any(w < 0))
    stop("negative weights not allowed")
  if (is.character(family))
    family <- get(family)
  if (is.function(family))
    family <- family()
  if (!inherits(family, "vglmff")) {
    stop("'family = ", family, "' is not a VGAM family function")
  }
  eval(vcontrol.expression)
  if (length(slot(family, "first")))
    eval(slot(family, "first"))
  vglm.fitter <- get(method)
  fit <- list(x = x, y = y, w = w, offset = offset,
              Xm2 = Xm2, Ym2 = Ym2, etastart = etastart, mustart = mustart,
              coefstart = coefstart, family = family, control = control,
              constraints = constraints, extra = extra, qr.arg = qr.arg,
              Terms = mt, function.name = function.name)
  return(fit)
}

