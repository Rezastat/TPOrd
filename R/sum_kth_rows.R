sum_kth_rows <- function(original_matrix, start_row, k) {
  if (!is.matrix(original_matrix)) {
    stop("The input original_matrix must be a matrix.")
  }
  if (start_row <= 0 || start_row > nrow(original_matrix)) {
    stop("The start_row is out of the matrix's bounds.")
  }
  if (k <= 0) {
    stop("Parameter k must be a positive integer.")
  }

  # keep rows before start_row
  if(start_row!=1){
    new_matrix <- original_matrix[1:(start_row - 1), ]
  } else{ new_matrix <- NULL}

  # sum every kth row starting from start_row
  for (i in seq(start_row, nrow(original_matrix), by = k)) {
    row_end <- min(i + k - 1, nrow(original_matrix))
    sum_row <- colSums(original_matrix[i:row_end, ])
    new_matrix <- rbind(new_matrix, sum_row)
  }

  return(new_matrix)
}
