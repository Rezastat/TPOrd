obsIM_R_ord <- function (formula, miscov, auxvar, model_type, data, beta, p_gz, num_categories)
{

  terms_dat <- unique(c(all.vars(formula), attr(terms(miscov),
                                                "term.labels"), attr(terms(auxvar), "term.labels")))
  Y_ <- as.character(formula[[2]])
  G_ <- all.vars(miscov)
  Z_ <- all.vars(auxvar)
  p_gz <- p_gz[p_gz$q > 0, ]
  g.vals <- unique(p_gz[, G_])
  N <- NROW(data)



  dat_ <- cbind(id_ = rep(seq(1, N, by = 1), each = NROW(g.vals)),
                data[rep(seq_len(N), each = NROW(g.vals)), ], matrix(rep(t(g.vals),
                                                                         N), ncol = NCOL(g.vals), byrow = TRUE))
  names(dat_)[(NCOL(dat_) - NCOL(g.vals) + 1):ncol(dat_)] <- G_
  dat_ <- droplevels(dat_)

  wgs <- omega1_g_ord(formula, model_type, dat_, beta, p_gz, G_, Z_)
  dat_$wg_ <- wgs$wg
  dat_ <- dat_[!is.na(dat_$wg_), ]


  IM1 <- obsIMf1_ord(formula,dat_,beta,p_gz,G_,Z_,model_type=model_type, num_categories= num_categories)$I3/N
  IM23_id <- obsIMf1_ord(formula,dat_,beta,p_gz,G_,Z_,model_type=model_type, num_categories= num_categories, by.id = TRUE)
  IM2_id <- IM23_id$I2/N
  IM3_id <- IM23_id$I3/N
  npars <- NCOL(IM1)
  Ind <- rep(FALSE, npars - 1)
  Betas_ids <- which(grepl(paste0(paste0("^", G_), collapse = "|"),
                           colnames(model.matrix(formula, dat_[1, ]))))
  Ind[Betas_ids] <- TRUE
  Null <- all(beta[Ind] == 0)
  return(list(IM1 = IM1, IM2_id = IM2_id, IM3_id = IM3_id,
              npars = npars, Ind = Ind, Null = Null))
}
