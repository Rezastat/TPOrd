twoPhaseSPML_ord <- function(formula, miscov, auxvar, family=NULL, data0, data1, start.values = NULL,
                             verbose = FALSE, n_second, model_type, num_categories, N, osm = T,
                             iteration = NULL, em_tol = 1e-05, max_iter = 500, array_id = NULL) {



  terms_dat <- unique(c(all.vars(formula), attr(terms(miscov), "term.labels"), attr(terms(auxvar), "term.labels")))
  Y_ <- as.character(formula[[2]])
  G_ <- all.vars(miscov)
  Z_ <- all.vars(auxvar)

  terms_data <- if(all(G_ %in% all.vars(formula))) {
    c(all.vars(formula)[-which(all.vars(formula) %in% G_)], Z_)
  } else {
    c(all.vars(formula), Z_)
  }
  n1 <- nrow(data1)
  n0 <- nrow(data0)
  n <- n1 + n0
  nG_ <- length(G_)
  uniqZ <- unique(data1[, Z_])
  uniqZ <- uniqZ[order(uniqZ)]
  uniqG <- unique(data1[, G_])

  dat0 <- cbind(id_ = rep(seq(n1 + 1, n, by = 1), each = NROW(uniqG)),
                data0[rep(seq_len(n0), each = NROW(uniqG)),],
                matrix(rep(t(uniqG), n0), ncol = NCOL(uniqG), byrow = TRUE))
  names(dat0)[(ncol(dat0) - NCOL(uniqG) + 1):ncol(dat0)] <- G_

  qform <- as.formula(paste0("wg_ ~", paste(G_, collapse = "+"), "+", Z_))

  dat <<- rbind(cbind(id_ = 1:n1, data1[, terms_dat]), dat0[, c("id_", terms_dat)])

  uniqGZ <- unique(data1[, c(G_, Z_)])
  q0 <- cbind(uniqGZ, q = 1 / nrow(uniqGZ))
  Zuniq <- unique(c(data1[, Z_], data0[, Z_]))

  q <- cbind(matrix(rep(t(uniqG), length(Zuniq)), ncol = NCOL(uniqG), byrow = TRUE), Z = rep(Zuniq, each = NROW(uniqG)))
  colnames(q) <- c(G_, Z_)
  q <- merge(q, q0, all.x = T, sort = F)
  q[is.na(q$q), "q"] <- 0

  nq <- nrow(q)
  if(!is.null(start.values$q)) {
    q1 <- merge(q, start.values$q, by = c(G_, Z_), all.x = T, sort = F)
    q1 <- q1[, c(G_, Z_, "q.y")]
    names(q1)[3] <- "q"
    q <- q1
  }

  X <- model.matrix(formula, dat)

  if(is.null(start.values$betas)) {
    if(model_type == "Stereotype_Regression") {
      mod_osm <- OSM_Weighted(convertResponseToFactor(formula), data = dat)
      theta.o <- c(mod_osm$phi[2:(dim(mod_osm$phi)[1]-1), 1] ,
                   mod_osm$alpha[2:(dim(mod_osm$alpha)[1]), 1] ,
                   mod_osm$beta[, 1])
      # theta.o <- c(mod_osm$phi[2:(dim(mod_osm$phi)[1]-1), 1] + runif(length(mod_osm$phi)-2, -0.0001, 0.0001),
      #              mod_osm$alpha[2:(dim(mod_osm$alpha)[1]), 1] + runif(length(mod_osm$alpha)-1, -0.0001, 0.0001),
      #              mod_osm$beta[, 1] + runif(length(mod_osm$beta), -0.0001, 0.0001))
    } else {
      mod_vglm <- vglm(formula, data = dat, family = family)
      theta.o <- mod_vglm@coefficients
    }
  } else {
    theta.o <- start.values$betas
  }

  nparms <- length(theta.o)
  ng <- NROW(uniqG)
  n1ones <- rep(1, n1)

  maxiter <- max_iter
  tol <- em_tol
  iter <- 0
  theta <- theta.o
  theta.o <- theta + 1
  q.o <- q
  q.o[, "q"] <- q.o[, "q"] - 1

  if(verbose) cat("Initial iteration:", iter, "theta =", theta, "pGZ =", as.numeric(q[, "q"]), "\n")

  while(iter <= maxiter & (max(abs(theta - theta.o)) > tol | max(abs(q[, "q"] - q.o[, "q"])) > tol)) {

    # if (iteration == 1 & array_id == 1){
    #   cat(sprintf("Iteration: %d theta_tol: %f q_tol: %f\n", iter, max(abs(theta - theta.o)), max(abs(q[, "q"] - q.o[, "q"]))), file = log_em, append = TRUE)
    # }
    cat("Iteration:", iter, "| theta_tol:", max(abs(theta - theta.o)) , "| q_tol:", max(abs(q[, "q"] - q.o[, "q"])),"\n")

    theta.o <- theta
    q.o <- q
    iter <- iter + 1

    # E-step
    dat$wg_ <<- c(n1ones, fnc_wgts_ordinal(formula, family=family, dat0, theta, q, G_, Z_, uniqZ, num_categories, N_total = N, n2 = n_second, model_type))

    # M-step
    if (model_type != "Stereotype_Regression" | osm == F ){
      dat$wg_[dat$wg_ == 0] <<- 1.0e-8
    }





    if(model_type == "Stereotype_Regression") {
      if (osm == T) {
        fit <- OSM_Weighted(formula = convertResponseToFactor(formula), data = dat, weight = dat$wg_)
      } else {
        fit <- rrvglm(formula = minus_to_response(formula), family = multinomial(), data = dat, weights = dat$wg_)
      }

    } else {
      fit <- vglm(formula = formula, family = family, data = dat, weights = dat$wg_)
    }

    theta <- if(model_type == "Stereotype_Regression") {
      if (osm == T){
        c(fit$phi[2:(dim(fit$phi)[1]-1), 1],
          fit$alpha[2:(dim(fit$alpha)[1]), 1],
          fit$beta[, 1])
      } else {
        rrvglmfit_sum <- summary(fit)
        adjust_coefficients(rrvglmfit_sum, num_categories)[, "Estimate"]
      }

    } else {
      c(fit@coefficients)
    }


    q.it <- data.frame(xtabs(qform, data = dat) / n)
    ordq <- match(apply(q[, c(G_, Z_)], 1, paste, collapse = "-"), apply(q.it[, c(G_, Z_)], 1, paste, collapse = "-"))
    q[, "q"] <- q.it[ordq, "Freq"]

    if(verbose) cat("Iteration:", iter, "theta =", theta, "pGZ =", as.numeric(q[, "q"]), "\n")
  }

  if(verbose) cat("Total iterations:", iter, "theta_end =", theta, "pGZ_end =", as.numeric(q[, "q"]), "\n")

  #Create objects needed in the testing part
  qG <- aggregate(as.formula(paste0("q~",paste(G_,collapse="+"))),q,FUN=sum)
  qG$k_ <- 1:nrow(qG)
  dat$ord_ <- 1:nrow(dat)
  dat <- merge(dat,qG[,c(G_,"k_")],by=G_) #order is messed up here
  dat <- dat[order(dat$ord_),]
  names(theta) <- gsub("beta_", "", names(theta))

  if( Ho<-!all(G_ %in% names(theta)) ){ ## under the null (no G in formula)

    indx_thetas_H0 <- which(!G_ %in% names(theta))
    nG_ <- length(indx_thetas_H0)
    val <- c(theta, rep(0,nG_))
    if (model_type == "Stereotype_Regression"){
      id <- c(seq_along(theta), (2*num_categories-3)+seq(0.1,0.9,length.out=nG_))
    } else {
      id <- c(seq_along(theta), (num_categories-1)+seq(0.1,0.9,length.out=nG_))
    }

    theta1 <- val[order(id)]
    names(theta1)[which(names(theta1)=="")] <- G_[indx_thetas_H0]
    formula1 <- update(formula,paste0(".~",paste(G_[indx_thetas_H0],collapse="+"),"+ ."))


    if (model_type != "Stereotype_Regression"){
      vglmfitparams <-vglm_qr(formula1, family = family,  data = dat, weights = dat$wg_, etastart = NULL, mustart = NULL, coefstart = theta1)

      vglmfit1 <- do.call(vglm_fit_qr,vglmfitparams)
    }


    #################
    if (model_type != "Stereotype_Regression"){
      FIM<-FIM_ordinal(theta1,qG,formula1,dat,fit=vglmfit1,model_type=model_type, n_second = n_second, num_categories = num_categories)
    } else {
      FIM<-FIM_ordinal(theta1,qG,formula1,dat,model_type=model_type, n_second = n_second, num_categories = num_categories)
    }



    res_SobsIM <- Testing_EM_joint0_ord(theta1,qG,Betas_ids=NULL,G_[indx_thetas_H0],FIM)


  }else{ ## under the alternative (G in formula)

    if (model_type != "Stereotype_Regression"){
      vglmfit0 <- list()
      vglmfit0[["tfit"]][["qr"]] <-fit@R
      FIM<-FIM_ordinal(theta,qG,formula,dat,fit=vglmfit0,model_type=model_type, n_second = n_second, num_categories = num_categories)
    } else {
      FIM<-FIM_ordinal(theta,qG,formula,dat,model_type=model_type, n_second = n_second, num_categories = num_categories)
    }



    res_VobsIM <- Testing_EM_joint0_ord(theta,qG,NULL,G_,FIM)

  }

  ll <- loglik_ord(theta,q,formula,Y_,G_,Z_,dat,family,num_categories,model_type) ## log-likelihood

  if(Ho){
    return(list(theta=theta,qGZ=q,Sobs=res_SobsIM$S,qG=qG,iter=iter,ll=ll,FIM=FIM$FIM))
  }else{
    return(list(theta=theta,qGZ=q,var_theta=res_VobsIM$Var,Wobs=res_VobsIM$W,qG=qG,iter=iter,ll=ll,FIM=FIM$FIM))
  }
}
