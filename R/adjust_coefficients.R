adjust_coefficients <- function(mos_str_sum, num_categories) {
  # get the the coefficients matrix
  coef_matrix <- mos_str_sum@coef3
  row_names <- row.names(coef_matrix)
  # Reversing the score parameters
  if (num_categories>3){
    score_param <- coef_matrix[1:(num_categories-2),]
    score_parameters <- score_param[(num_categories-2):1,]
  } else {
    score_parameters <- coef_matrix[1, ]
  }
  intercepts <- coef_matrix[(2 * num_categories-3):(num_categories-2+1),]

  coefs <- coef_matrix[(2 * num_categories-2):nrow(coef_matrix),]

  final_matrix <- rbind(score_parameters,intercepts,coefs)
  row.names(final_matrix) <- row_names

  return(final_matrix)
}
