generateMaxCorrelatedData <- function(n, fractions_cat1, fractions_cat2, negative = FALSE, rand = FALSE) {
  if (sum(fractions_cat1) > 1 || sum(fractions_cat2) > 1) {
    stop("The sum of fractions must not exceed 1.")
  }
  # random data generation
  if (rand) {
    # Generating random data based on fractions
    variable1 <- sample(1:(length(fractions_cat1) + 1), size = n, prob = c(fractions_cat1, 1 - sum(fractions_cat1)), replace = T)
    variable2 <- sample(1:(length(fractions_cat2) + 1), size = n, prob = c(fractions_cat2, 1 - sum(fractions_cat2)), replace = T)


  } else {
    # deterministic data generation
    variable1 <- rep(NA, n)
    variable2 <- rep(NA, n)

    ends_cat1 <- cumsum(c(fractions_cat1, 1 - sum(fractions_cat1))) * n
    ends_cat2 <- cumsum(c(fractions_cat2, 1 - sum(fractions_cat2))) * n

    for (i in 1:(length(ends_cat1) - 1)) {
      variable1[(ifelse(i == 1, 1, ends_cat1[i-1] + 1)):ends_cat1[i]] <- i
    }
    variable1[(ends_cat1[length(ends_cat1)-1] + 1):n] <- length(ends_cat1)

    for (j in 1:(length(ends_cat2) - 1)) {
      variable2[(ifelse(j == 1, 1, ends_cat2[j-1] + 1)):ends_cat2[j]] <- j
    }
    variable2[(ends_cat2[length(ends_cat2)-1] + 1):n] <- length(ends_cat2)
  }

  # Sorting for max correlation
  variable1 <- sort(variable1)
  variable2 <- if (negative) { sort(variable2, decreasing = TRUE) } else { sort(variable2) }

  max_cor <- cor(as.numeric(variable1), as.numeric(variable2))

  # cat(sprintf("Correlation achieved: %f\n", max_cor))

  return(list(variable1 = variable1, variable2 = variable2, correlation = max_cor))
}



