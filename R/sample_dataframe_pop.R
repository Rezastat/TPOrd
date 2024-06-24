sample_dataframe_pop <- function(dat, n2, m, j) {
  # Define unique categories
  num_categories_y <- length(unique(dat$Y))
  num_categories_z <- length(unique(dat$Z))
  dat$YZ <- interaction(dat$Y, dat$Z)
  num_categories_yz <- length(unique(dat$YZ))

  dat$R1 <- 0
  dat$R2 <- 0
  dat$R3 <- 0
  dat$R4 <- 0

  n_fixed <- floor(n2 * m)

  # R1: Simple Random Sampling
  random_sample_indices <- sample(nrow(dat), n_fixed)
  dat$R1[random_sample_indices] <- 1

  # R2: Outcome-dependent sampling
  dat <- TPOrd:::apply_balanced_sampling(dat, "Y", n_fixed, num_categories_y, "R2")
  R2_indices <- which(dat$R2 ==1)
  # R3: Covariate-dependent sampling
  dat <- TPOrd:::apply_balanced_sampling(dat, "Z", n_fixed, num_categories_z, "R3")
  R3_indices <- which(dat$R3 ==1)
  # R4: Covariate and outcome-dependent sampling
  dat <- TPOrd:::apply_balanced_sampling(dat, "YZ", n_fixed, num_categories_yz, "R4")
  R4_indices <- which(dat$R4 ==1)

  # Create empty matrices for R2, R3, and R4
  R2_matrix <- matrix(0, nrow = nrow(dat), ncol = j)
  R3_matrix <- matrix(0, nrow = nrow(dat), ncol = j)
  R4_matrix <- matrix(0, nrow = nrow(dat), ncol = j)

  # Perform balanced sampling for R2, R3, and R4 while keeping the fixed sample
  for (iteration in 1:j) {
    R2_matrix[, iteration] <- apply_balanced_sampling_fixed(dat, "Y", n2, num_categories_y, "R2", R2_indices)$R2
    R3_matrix[, iteration] <- apply_balanced_sampling_fixed(dat, "Z", n2, num_categories_z, "R3", R3_indices)$R3
    R4_matrix[, iteration] <- apply_balanced_sampling_fixed(dat, "YZ", n2, num_categories_yz, "R4", R4_indices)$R4
  }

  return(list(dat = dat, R2_matrix = R2_matrix, R3_matrix = R3_matrix, R4_matrix = R4_matrix))
}

apply_balanced_sampling_fixed <- function(dat, strat_var, n2, num_categories, r_var, fixed_indices) {
  category_counts <- table(dat[[strat_var]])
  categories <- names(category_counts)

  samples_per_category <- rep(0, length(categories))
  names(samples_per_category) <- categories

  # Initial allocation of samples ensuring not to exceed the per-category count
  for (cat in categories) {
    proposed_allocation <- floor(n2 / num_categories)
    samples_per_category[cat] <- min(proposed_allocation, category_counts[cat])
  }

  remaining_samples <- n2 - sum(samples_per_category)
  remaining_categories <- categories[samples_per_category < category_counts]

  # Redistribute remaining samples
  while (remaining_samples > 0) {
    remaining_categories <- categories[samples_per_category < category_counts]
    if (length(remaining_categories) <= remaining_samples) {
      # If there are fewer or equal remaining categories than samples, add one to each
      samples_per_category[remaining_categories] <- samples_per_category[remaining_categories] + 1
      remaining_samples <- remaining_samples - length(remaining_categories)
    } else {
      # More remaining samples than categories, randomly select categories to increment
      categories_to_increment <- sample(remaining_categories, remaining_samples, replace = FALSE)
      samples_per_category[categories_to_increment] <- samples_per_category[categories_to_increment] + 1
      remaining_samples <- 0
    }
  }

  # Mark all fixed indices
  dat[[r_var]] <- 0
  dat[[r_var]][fixed_indices] <- 1

  # Exclude fixed indices from available sampling pool
  available_indices <- setdiff(seq_len(nrow(dat)), fixed_indices)

  for (cat in categories) {
    category_indices <- which(dat[[strat_var]] == cat)
    available_category_indices <- intersect(category_indices, available_indices)
    if (length(available_category_indices) == 1) {
      available_category_indices <- c(available_category_indices, available_category_indices)
    }

    num_samples_to_assign <- samples_per_category[cat] - sum(dat[[r_var]][category_indices])

    if (num_samples_to_assign > 0) {
      if (length(available_category_indices) >= num_samples_to_assign) {
        # Perform sampling
        stratified_sample_indices <- sample(available_category_indices, num_samples_to_assign, replace = FALSE)
        dat[[r_var]][stratified_sample_indices] <- 1
      } else {
        stop(paste("Error: Not enough indices available in category", cat,
                   "to assign", num_samples_to_assign, "samples.",
                   "Available indices:", length(available_category_indices)))
      }
    }
  }

  return(dat)
}





