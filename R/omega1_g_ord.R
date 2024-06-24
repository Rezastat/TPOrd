omega1_g_ord <- function(formula, model_type, dat0, beta, q0, G_, Z_) {
  N <- wg0 <- id_ <- lwg0 <- NULL # Initialize to avoid global binding issues

  resp <- as.character(formula[[2]])
  terms_q <- c(G_, Z_)
  dat0$ord <- 1:nrow(dat0)

  q_ <- names(q0)[!names(q0) %in% terms_q]
  if (q_ %in% names(dat0)) {
    namsold <- names(dat0)[which(names(dat0) == q_)]
    names(dat0)[which(names(dat0) == q_)] <- paste0(q_, ".old")
  }

  dat0bis <- merge(as.data.table(dat0), as.data.table(q0[, c(terms_q, "q")]), by = terms_q, all.x = T, sort = F)

  X0 <- model.matrix(formula, dat0bis)
  extended_X0 <- expand_model_matrix(X0, num_categories - 1)

  big_N <- nrow(dat0)

  if (model_type == "Stereotype_Regression") {
    phi_col <- rep(c(beta[1:(num_categories - 2)], 1), big_N)
    extended_X0[, num_categories:ncol(extended_X0)] <- extended_X0[, num_categories:ncol(extended_X0)] * phi_col
    eta_vals <- extended_X0 %*% beta[(num_categories - 1):length(beta)]
  } else {
    eta_vals <- extended_X0 %*% beta
  }

  eta_matrix <- matrix(eta_vals, nrow = big_N, ncol = num_categories - 1, byrow = T)

  if (model_type == "Proportional_Odds") {
    cum_probs <- exp(eta_matrix) / (1 + exp(eta_matrix))
    p_vals <- matrix(nrow = big_N, ncol = num_categories)
    p_vals[, 1] <- cum_probs[, 1]

    for (j in 2:(num_categories - 1)) {
      p_vals[, j] <- cum_probs[, j] - cum_probs[, j - 1]
    }

    p_vals[, num_categories] <- 1 - cum_probs[, num_categories - 1]

  } else if (model_type == "Adjacent_Category") {
    exp_cum_sum_matrix <- exp(t(apply(eta_matrix, 1, function(x) rev(cumsum(rev(x))))))
    p_vals <- matrix(nrow = big_N, ncol = num_categories)

    for (i in 1:big_N) {
      denominator <- 1 + sum(exp_cum_sum_matrix[i, ])
      for (j in 1:(num_categories - 1)) {
        p_vals[i, j] <- exp_cum_sum_matrix[i, j] / denominator
      }
      p_vals[i, num_categories] <- 1 / denominator
    }

  } else if (model_type == "Stopping_Ratio") {
    exp_sum_matrix <- exp(eta_matrix)
    p_vals <- matrix(nrow = big_N, ncol = num_categories)

    for (i in 1:big_N) {
      for (j in 1:(num_categories - 1)) {
        if (j == 1) {
          p_vals[i, j] <- exp_sum_matrix[i, j] / (1 + exp_sum_matrix[i, j])
        } else {
          product_term <- 1
          for (r in 1:(j - 1)) {
            product_term <- product_term * (1 - exp_sum_matrix[i, r] / (1 + exp_sum_matrix[i, r]))
          }
          p_vals[i, j] <- exp_sum_matrix[i, j] / (1 + exp_sum_matrix[i, j]) * product_term
        }
      }
    }
    p_vals[, num_categories] <- 1 - rowSums(p_vals[, 1:(num_categories - 1)])

  } else if (model_type == "Stereotype_Regression") {
    exp_sum_matrix <- cbind(rep(1, big_N), exp(eta_matrix))
    p_vals <- matrix(nrow = big_N, ncol = num_categories)
    for (j in 1:num_categories) {
      p_vals[, j] <- exp_sum_matrix[, j] / rowSums(exp_sum_matrix)
    }
  }

  inner_arg <- sapply(1:big_N, function(i) {
    probs_i <- p_vals[i, ]
    x_vec <- integer(num_categories)
    x_vec[dat0bis$Y[i]] <- 1
    dmultinom(x = x_vec, size = 1, prob = probs_i)
  })

  dat0bis$lwg0 <- log(inner_arg) + log(dat0bis[["q"]])
  dat0bis$wg0 <- exp(dat0bis$lwg0)
  dat0bis[, N := sum(wg0, na.rm = TRUE), by = id_] # Sum of weights by id

  if (any(dat0bis$N == 0)) {
    zerindx <- which(dat0bis$N == 0)
    dfwg0 <- dat0bis[zerindx, c("id_", "ord", "lwg0"), with = FALSE]
    dfwg0[, max := max(lwg0), by = id_]
    if (!all(dfwg0$ord == dat0bis$ord[zerindx])) stop("Something went wrong in the probability normalization process.")
    dat0bis$wg0[zerindx] <- exp(dfwg0$lwg0 - dfwg0$max)
    dat0bis[zerindx, N := sum(wg0), by = id_]
  }

  dat0bis$wg <- dat0bis$wg0 / dat0bis$N
  dat0bis <- dat0bis[order(dat0bis$ord),]

  return(list(wg = dat0bis$wg, pYZ = dat0bis$N))
}
