expand_model_matrix <- function(model_matrix, num_cat) {
  if (num_cat < 2) {
    stop("num_cat should be greater or equal to 2")
  }

  # intercept columns
  intercept_cols <- sapply(1:num_cat, function(i) {
    col <- rep(0, nrow(model_matrix) * num_cat)
    col[seq(i, length(col), num_cat)] <- 1
    col
  })
  colnames(intercept_cols) <- paste0("(Intercept):", 1:num_cat)

  # Replicate rows
  replicated_matrix <- model_matrix[rep(1:nrow(model_matrix), each = num_cat), ]

  # Combine intercept with the replicated matrix
  expanded_matrix <- cbind(intercept_cols, replicated_matrix[, -1])

  # row names
  rownames(expanded_matrix) <- paste0(rep(rownames(model_matrix), each = num_cat), ":", rep(1:num_cat, nrow(model_matrix)))

  return(expanded_matrix)
}
