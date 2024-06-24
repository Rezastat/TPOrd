change_and_copy_names <- function(original_vector, target_vector) {

  if (length(target_vector) != length(original_vector)) {
    stop("The target vector must be of the same length as the original vector.")
  }

  # Get the names of the original vector
  original_names <- names(original_vector)

  # Replace 'est' with 'se' in the names
  new_names <- gsub("est", "se", original_names)

  # Assign the new names to the target vector
  names(target_vector) <- new_names


  return(target_vector)
}
