capture_fit_with_warnings <- function(Y, Z, model_type) {
  warningOccurred <- FALSE
  lmfit.comp <- tryCatch({
    withCallingHandlers({
      if (model_type == "Proportional_Odds") {
        vglm(Y ~ Z, family = cumulative(parallel = TRUE, reverse = FALSE))
      } else if (model_type == "Adjacent_Category") {
        vglm(Y ~ Z, family = acat(reverse = TRUE, parallel = TRUE))
      } else if (model_type == "Stopping_Ratio") {
        vglm(Y ~ Z, family = sratio(reverse = FALSE, parallel = TRUE))
      } else if (model_type == "Stereotype_Regression") {
        tempdf <- data.frame(Y = factor(Y), Z = Z)

        OSM_Weighted(Y ~ Z, data = tempdf, Hess = TRUE)

      } else {
        stop("Invalid model type")
      }
    }, warning = function(w) {
      if (grepl("NaNs produced", w$message)) {
        warningOccurred <<- TRUE
        message("A warning of interest occurred: ", w$message)
      }
    })
  }, error = function(e) {
    message("An error occurred during model fitting: ", e$message)
    return(NULL)
  })

  list(lmfit.comp = lmfit.comp, warningOccurred = warningOccurred)
}
