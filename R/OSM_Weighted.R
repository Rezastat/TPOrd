# Program to estimate the parameters of the Stereotype Model and Ordered Stereotype Model using command "optim"
# Author: Laia Egea
# Program created: 04-05-2022
# Last modification: 15-11-2022

# library(msm) # deltamethod

# p = number of covariates
# q = number of categories resp. var.
# phi1.iszero = TRUE means phi_1=0 ; FALSE means phi_1=1

########################## Functions needed in the final function ##########################

# The inverse of the logit function
expit <- function(x){
  1/(1+exp(-x))
}


# Function to define the vector of estimated parameters in a list:
unpack.par.stereo <- function(param, p, q, phi1.iszero){

  alpha <- c(0,param[1:(q-1)])

  if (phi1.iszero == TRUE) {
    phi <- c(0, param[(q-1)+1:(q-2)], 1)
  } else {
    phi <- c(1, param[(q-1)+1:(q-2)], 0)
  }

  beta <- param[q+q-3+1:p]

  return(list(alpha=alpha, phi=phi, beta=beta))
}


# Function to compute the log(L) for each individual (i) and each level of the response variable (k)

compute.elem.loglike.stereo <- function(k,i,matrixdata, parlist, weight){

  NumCol <- ncol(matrixdata)

  #First, the product between effects and covariates (beta*x_i)
  effectcov <- sum(parlist$beta * matrixdata[i,2:NumCol])

  #Second, the elements of the loglikelihood
  element1 <- parlist$alpha[k] + (parlist$phi[k]*effectcov)
  element2 <- sum(exp(parlist$alpha+(parlist$phi*effectcov)))
  element <- (element1 - log(element2))

  return(element)
}

# Function to reverse D.Fernandez reparametrization (Ordered Stereotype Model)
inverse.repar.param <- function(parlist, phi1.iszero, q){

  if(q == 3){
    parlist$phi[2] <- expit(parlist$phi[2])
  } else{
    if(phi1.iszero){
      parlist$phi[2:(q-1)] <- expit(cumsum(c(parlist$phi[2], exp(parlist$phi[3:(q-1)]))))
    } else {
      parlist$phi[2:(q-1)] <- expit(cumsum(c(parlist$phi[2], -exp(parlist$phi[3:(q-1)]))))
    }
  }

  return(parlist$phi)
}

# Set up the function call to use in optim(), which extracts the parameters from the parameter vector
loglikelihood.stereo <- function(param, matrixdata, lpar, phi1.iszero, repar, weight = weight)
{
  # Define some parameters
  p <- lpar$p #number of explanatory variables
  q <- lpar$q #number of category data
  N <- nrow(matrixdata)

  # Transform de vector of parameter estimates to list
  parlist <- unpack.par.stereo(param, p, q, phi1.iszero)

  # Transform the reparametrization
  if(repar){
    parlist$phi <- inverse.repar.param(parlist, phi1.iszero, q)
  }

  # Compute the log(L) for each individual (i) and his level of the response variable (k)
  if (length(weight)!=0){
    loglike.stereo <- sum(weight * sapply(1:N, function(i) compute.elem.loglike.stereo(k = matrixdata[i,"Y"],
                                                                                       i = i,
                                                                                       matrixdata = matrixdata,
                                                                                       parlist = parlist)))
  } else{
    loglike.stereo <- sum(sapply(1:N, function(i) compute.elem.loglike.stereo(k = matrixdata[i,"Y"],
                                                                              i = i,
                                                                              matrixdata = matrixdata,
                                                                              parlist = parlist)))
  }

  return(-loglike.stereo)
}



########################## Function to fit the SM and OSM ##########################

OSM_Weighted <- function(formula, data,  phi1.iszero = T, repar = T, Approx.start.values = FALSE, method = "BFGS",
                         control = list(trace=TRUE, REPORT=1, maxit=10000, reltol=1e-15), weight = NULL, Hess = FALSE, ...){

  response <- all.vars(formula[[2]])
  covariates <- all.vars(formula[[3]])

  # Transform the data to estimate the model: create the dummy variable if it is needed
  matrixdata <- cbind(Y = data[,response], model.matrix(as.formula(paste0("~ ",paste(covariates, collapse = " + "))), data))
  matrixdata <- matrixdata[,-which(colnames(matrixdata) == "(Intercept)")]

  if(sum(apply(matrixdata,2, function(x) any(is.na(x))))>0){
    matrixdata <- matrixdata[complete.cases(matrixdata),]
  }

  if(min(matrixdata[,1]) == 0){
    matrixdata[,1] <- matrixdata[,1]+1
  }


  # Keep the real names of the covariates
  OUTnames <- colnames(matrixdata)[-1]

  # Create a list with some important parameters:
  lpar <- list()
  lpar$n <- nrow(matrixdata) # number of individuals
  lpar$p <- ncol(matrixdata) - 1 # number of covariates
  lpar$q <- length(table(matrixdata[,"Y"])) # number of categories resp. var.
  lpar$numParam <- (2*lpar$q) - 3 + lpar$p # number of parameters to estimate

  # Trnasform the names: Y = response variable and Xs = for covariates
  colnames(matrixdata) <- c("Y", paste0("X", 1:lpar$p))


  # Rename some parameters
  p <- lpar$p
  q <- lpar$q
  numParam <- lpar$numParam

  # Initialize parameters
  if(Approx.start.values){
    n <- lpar$n

    ############# BETAS #############
    # try logistic/probit regression on 'middle' cut to find starting
    # values for the coefficients of the predictors
    # q1 is the level at, or just before, halfway through the levels of y
    (q1 <- q %/% 2L)

    ## y1 is a binary variable with y1=0 if y <= q1 and y1=1 if y > q1
    y1 <- (matrixdata[,1] > q1)
    table(y1)

    ## Construct a new model matrix and add an intercept column to it
    X <- cbind(Intercept = rep(1, n), matrixdata[,-1])
    head(X)

    ## Now attempt to fit logistic regression to the binary response y1
    fit <- glm.fit(X, y1, family = binomial())
    if(!fit$converged)
      stop("attempt to find suitable starting values failed")
    coefs <- fit$coefficients
    if(any(is.na(coefs))) {
      warning("design appears to be rank-deficient, so dropping some coefs")
      keep <- names(coefs)[!is.na(coefs)]
      coefs <- coefs[keep]
      x <- x[, keep[-1L], drop = FALSE]
      pc <- ncol(x)
    }

    ############# ALPHAS #############
    ## The other parameters are labelled as alphas in Agresti's definition
    ## of the proportional odds model. They are the base probabilities
    ## for each level of the response variable, and must be strictly increasing.
    ## Generate them initially assuming they're evenly spaced across the
    ## range, convert to the linear predictor space using the logit link
    ## and adjust them to incorporate the fact that the logistic
    ## regression fitting produced an intercept term coefs[1L] for the
    ## q1 level
    logit <- function(p) log(p/(1 - p))
    spacing <- logit((1L:(q-1))/((q-1)+1L)) # just a guess
    gammas <- -coefs[1L] + spacing - spacing[q1]

    ############# PHIs #############
    ## Also generate starting values for phi, assuming equal spacing in
    ## the space of phi and converting using the logit link to the space
    ## of the auxiliary variable u
    startingphi <- (1:(q-2))/(q-1)
    u2 <- logit(startingphi)[1]
    us <- log(diff(logit(startingphi)))


    ############# JOINT #############
    ## Construct the full starting values vector, using the fact that
    ## coefs[1L] has already been incorporated into the gammas object
    par.start <- c(gammas, u2, us, coefs[-1L])
  } else {
    par.start <- rep(1e-1, numParam)
  }

  # Run optim
  est.R <- optim(par = par.start, # Initialize parameters
                 matrixdata = matrixdata, # Data matrix
                 lpar = lpar, # Parameters
                 phi1.iszero = phi1.iszero, # TRUE means phi_1=0 ; FALSE means phi_1=1
                 repar = repar, # FALSE means nominal data, TRUE means ordinal data (reparametrization of the scores)
                 fn = loglikelihood.stereo, # Function to optimize
                 method = method, # Optimization method
                 hessian = Hess,
                 weight = weight,
                 ...) ### que devuelva la hessiana

  # Extract the results
  param.est.R <- unpack.par.stereo(est.R$par, p, q, phi1.iszero)



  # SE
  if (Hess) {
    VarCov <- solve(est.R$hessian)
    SE.R <- sqrt(diag(VarCov))
    param.SE.R <- unpack.par.stereo(SE.R, p, q, phi1.iszero)
  }


  # Reverse the reparametrization and apply the Delta Method to compute the SE of the scores
  if (repar){

    param.est.R$phi_norepar <- param.est.R$phi
    if(Hess) {
      param.SE.R$SEphi_norepar <- param.SE.R$phi
    }


    param.est.R$phi <- inverse.repar.param(param.est.R,phi1.iszero,q)

    # Delta method
    if (Hess) {
      if(q == 3){
        param.SE.R$phi[2] <- deltamethod(~ 1/(1+exp(-x1)), param.est.R$phi[2] , param.SE.R$phi[2])
      } else {
        Ind <- 2:(q-2)
        Ind <- paste0("x", Ind)

        if(phi1.iszero){
          Ind <- paste0("~1/(1 + exp(-((1/(1+exp(-x1))) + exp(", paste(Ind, collapse = ") + exp("), "))))")
          if(q >= 5){
            for(i in (q-3):2){
              Ind_aux <- 2:i
              Ind_aux <- paste0("x", Ind_aux)
              Ind_aux <- paste0("~1/(1 + exp(-((1/(1+exp(-x1))) + exp(", paste(Ind_aux, collapse = ") + exp("), "))))")
              Ind <- c(Ind_aux, Ind)
            }
            for(i in 2:(q-3)){
              for(j in 2:i){
                Ind[i] <- gsub(paste0("x",j),gsub("~", "", Ind[j-1]), Ind[i])
              }
            }
          }
        } else {
          Ind <- paste0("~1/(1 + exp(-((1/(1+exp(-x1))) - exp(", paste(Ind, collapse = ") - exp("), "))))")
          if(q >= 5){
            for(i in (q-3):2){
              Ind_aux <- 2:i
              Ind_aux <- paste0("x", Ind_aux)
              Ind_aux <- paste0("~1/(1 + exp(-((1/(1+exp(-x1))) - exp(", paste(Ind_aux, collapse = ") - exp("), "))))")
              Ind <- c(Ind_aux, Ind)
            }
            for(i in 2:(q-3)){
              for(j in 2:i){
                Ind[i] <- gsub(paste0("x",j),gsub("~", "", Ind[j-1]), Ind[i])
              }
            }
          }
        }

        Ind <- c("~ 1/(1+exp(-x1))", Ind)
        L <- lapply(Ind, as.formula)


        VarCov <- VarCov[(q-1)+1:(q-2), (q-1)+1:(q-2)]
        diag(VarCov) <- sqrt(diag(VarCov))

        param.SE.R$phi[2:(q-1)] <- deltamethod(L, param.est.R$phi[2:(q-1)] , VarCov)
      }
    }

  }


  deviance <- 2 * est.R$value

  results.R <- list()

  if (Hess) {
    results.R$alpha <- cbind(param.est.R$alpha, param.SE.R$alpha)
    colnames(results.R$alpha) <- c("est", "SE")
    rownames(results.R$alpha) <- paste0("alpha",1:q)
  } else {
    results.R$alpha <- matrix(param.est.R$alpha)
    colnames(results.R$alpha) <- "est"
    rownames(results.R$alpha) <- paste0("alpha",1:q)
  }


  if (Hess) {
    results.R$phi <- cbind(param.est.R$phi, param.SE.R$phi)
    results.R$phi[q,2] <- 0
    colnames(results.R$phi) <- c("est","SE")
    rownames(results.R$phi) <- paste0("phi",1:q)
  } else {
    results.R$phi <- matrix(param.est.R$phi)
    colnames(results.R$phi) <- "est"
    rownames(results.R$phi) <- paste0("phi",1:q)
  }

  if (Hess) {
    results.R$beta <- cbind(param.est.R$beta, param.SE.R$beta)
    colnames(results.R$beta) <- c("est","SE")
    rownames(results.R$beta) <- paste0("beta_",OUTnames)
  } else {
    results.R$beta <- matrix(param.est.R$beta)
    colnames(results.R$beta) <- "est"
    rownames(results.R$beta) <- paste0("beta_",OUTnames)
  }


  results.R$logLike <- est.R$value

  results.R$AIC <- 2*numParam + 2*results.R$logLike
  results.R$BIC <- 2*results.R$logLike + numParam*log(lpar$n)

  return(results.R)
}
