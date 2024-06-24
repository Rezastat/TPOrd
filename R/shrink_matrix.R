shrink_matrix <- function(original_matrix,num_categories){
  num_rows_new_matrix <- nrow(original_matrix) / (num_categories - 1)


  new_matrix <- matrix(0, nrow = num_rows_new_matrix, ncol = ncol(original_matrix))


  for (i in 1:num_rows_new_matrix) {
    row_start <- (i - 1) * (num_categories - 1) + 1
    row_end <- i * (num_categories - 1)
    new_matrix[i, ] <- colSums(original_matrix[row_start:row_end, ])
  }

  # keep column names
  colnames(new_matrix) <- colnames(original_matrix)
  return(new_matrix)
}
