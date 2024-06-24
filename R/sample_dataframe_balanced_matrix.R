sample_dataframe_balanced_matrix <- function(dat, n2, j) {
  # Define unique categories
  num_categories_y <- length(unique(dat$Y))
  num_categories_z <- length(unique(dat$Z))
  dat$YZ <- interaction(dat$Y, dat$Z)
  num_categories_yz <- length(unique(dat$YZ))

  # Initialize matrices to hold multiple versions of each sampling scheme
  R2_matrix <- matrix(0, nrow = nrow(dat), ncol = j)
  R3_matrix <- matrix(0, nrow = nrow(dat), ncol = j)
  R4_matrix <- matrix(0, nrow = nrow(dat), ncol = j)

  # Populate the matrices with samples
  for (i in 1:j) {
    # R2: Outcome-dependent sampling
    R2_matrix[, i] <- apply_balanced_sampling_matrix(dat, "Y", n2, num_categories_y)

    # R3: Covariate and outcome-dependent sampling
    R3_matrix[, i] <- apply_balanced_sampling_matrix(dat, "YZ", n2, num_categories_yz)

    # R4: Covariate-dependent sampling
    R4_matrix[, i] <- apply_balanced_sampling_matrix(dat, "Z", n2, num_categories_z)


  }

  return(list(R2 = R2_matrix, R3 = R3_matrix, R4 = R4_matrix))
}

apply_balanced_sampling_matrix <- function(dat, strat_var, n2, num_categories) {
  category_counts <- table(dat[[strat_var]])
  categories <- names(category_counts)

  samples_per_category <- rep(0, length(categories))
  names(samples_per_category) <- categories

  # Initial allocation of samples, ensuring not to exceed the per-category count
  for (cat in categories) {
    proposed_allocation <- floor(n2 / num_categories)
    samples_per_category[cat] <- min(proposed_allocation, category_counts[cat])
  }

  remaining_samples <- n2 - sum(samples_per_category)

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

  sample_vector <- rep(0, nrow(dat))

  for (cat in categories) {
    category_indices <- which(dat[[strat_var]] == cat)
    num_samples_to_assign <- samples_per_category[cat]

    if (num_samples_to_assign > 0) {
      if (num_samples_to_assign == 1){
        sample_vector[category_indices] <- 1
      } else{

        # Perform sampling
        stratified_sample_indices <- sample(category_indices, num_samples_to_assign, replace = FALSE)
        sample_vector[stratified_sample_indices] <- 1
      }
    }
  }

  return(sample_vector)
}
