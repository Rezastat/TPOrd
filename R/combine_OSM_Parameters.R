combine_OSM_Parameters <- function(osmList, num_categories) {
  # Extracting rows 2 to num_categories-1 from $phi
  phi_part <- osmList$phi[2:(num_categories-1), , drop = FALSE]
  phi_names <- attr(phi_part, "dimnames")[[1]]

  # Extracting rows 2 to end from $alpha
  alpha_part <- osmList$alpha[2:nrow(osmList$alpha), , drop = FALSE]
  alpha_names <- attr(alpha_part, "dimnames")[[1]]

  # Extracting all rows from $beta and its names
  beta_part <- osmList$beta
  beta_names <- attr(beta_part, "dimnames")[[1]]

  # Removing "beta_" from the beta names
  beta_names <- gsub("beta_", "", beta_names)

  # Combining the parts into one vector
  combined_vector <- c(phi_part, alpha_part, beta_part)

  # Setting names for the combined vector
  combined_names <- c(phi_names, alpha_names, beta_names)
  names(combined_vector) <- combined_names

  return(combined_vector)
}

