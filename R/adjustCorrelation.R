adjustCorrelation <- function(variable1, variable2, target_cor) {
  if (abs(target_cor) > 1) {
    stop("Target correlation must be between -1 and 1.")
  }

  current_cor <- cor(as.numeric(variable1), as.numeric(variable2))

  if (abs(current_cor) < abs(target_cor)) {
    stop("Target correlation is higher than the current maximum achievable correlation. Cannot increase correlation further.")
  }

  if (sign(current_cor) != sign(target_cor)) {
    stop("Current correlation and target correlation must have the same sign.")
  }

  iterations <- 0
  tol = 0.001
  # Determine if we are increasing or decreasing towards the target
  increasing <- target_cor > current_cor

  # Loop to reach the target for the specified direction
  while((increasing && current_cor < target_cor) || (!increasing && current_cor > target_cor)) {
    iterations <- iterations + 1

    idx1 <- sample(1:length(variable2), 1)
    idx2 <- sample(1:length(variable2), 1)

    temp <- variable2[idx1]
    variable2[idx1] <- variable2[idx2]
    variable2[idx2] <- temp

    new_cor <- cor(as.numeric(variable1), as.numeric(variable2))

    # Accept the swap if it's within the tolerance level of the target, even if slightly overshooting
    if ((increasing && new_cor > current_cor && new_cor <= target_cor + tol) ||
        (!increasing && new_cor < current_cor && new_cor >= target_cor - tol)) {
      current_cor <- new_cor  # Accept the swap
    } else {
      # Swap back: it didn't help or passed the tolerance
      variable2[idx2] <- variable2[idx1]
      variable2[idx1] <- temp
    }

    # Break condition based on absolute difference
    if (abs(current_cor - target_cor) <= tol) {
      break
    }

    if(iterations > 100000) {
      cat("Maximum iterations reached without achieving target correlation.\n")
      break
    }
  }

  # cat(sprintf("Target correlation approximately achieved: %f after %d iterations.\n", current_cor, iterations))

  return(list(variable1 = variable1, variable2 = variable2, correlation = current_cor))
}

