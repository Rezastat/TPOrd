sample_dataframe <- function(dat, n2) {
  # Define unique categories
  num_categories_y <- length(unique(dat$Y))
  num_categories_z <- length(unique(dat$Z))
  dat$YZ <- interaction(dat$Y, dat$Z)
  num_categories_yz <- length(unique(dat$YZ))

  dat$R1 <- 0
  dat$R2 <- 0
  dat$R3 <- 0
  dat$R4 <- 0

  # R1: Simple Random Sampling
  random_sample_indices <- sample(nrow(dat), n2)
  dat$R1[random_sample_indices] <- 1

  # R2: Outcome-dependent sampling
  dat <- apply_balanced_sampling(dat, "Y", n2, num_categories_y, "R2")

  # R3: Covariate-dependent sampling
  dat <- apply_balanced_sampling(dat, "Z", n2, num_categories_z, "R3")

  # R4: Covariate and outcome-dependent sampling
  dat <- apply_balanced_sampling(dat, "YZ", n2, num_categories_yz, "R4")

  return(dat)
}

apply_balanced_sampling <- function(dat, strat_var, n2, num_categories, r_var) {
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


  for (cat in categories) {
    category_indices <- which(dat[[strat_var]] == cat)
    if (length(category_indices) == 1){
      category_indices <- c(category_indices, category_indices)
    }


    num_samples_to_assign = samples_per_category[cat]

    if (num_samples_to_assign > 0) {
      if (length(category_indices) >= num_samples_to_assign) {
        # Perform sampling
        stratified_sample_indices = sample(category_indices, num_samples_to_assign, replace = FALSE)

        dat[[r_var]][stratified_sample_indices] = 1
      } else {
        stop(paste("Error: Not enough indices available in category", cat,
                   "to assign", num_samples_to_assign, "samples.",
                   "Available indices:", length(category_indices)))
      }
    }

  }



  return(dat)
}
