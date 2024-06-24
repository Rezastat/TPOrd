minus_to_response <- function(formula) {

  formula_parts <- as.character(formula)


  if (length(formula_parts) < 3 || formula_parts[1] != "~") {
    stop("Invalid formula format. Expected format: 'response ~ predictors'")
  }


  formula_parts[2] <- paste0("-", formula_parts[2])


  modified_formula_str <- paste(formula_parts[2], formula_parts[1], formula_parts[3])


  as.formula(modified_formula_str)
}
