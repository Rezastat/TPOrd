convertResponseToFactor <- function(formula) {
  response <- all.vars(formula)[1]
  predictors <- formula[[3]]

  # Create a new formula with the response variable as a factor
  newFormula <- as.formula(paste0("factor(", response, ") ~ ", deparse(predictors)))

  return(newFormula)
}
