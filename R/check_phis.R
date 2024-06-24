check_phis <- function(vec) {
  # see if all elements are between 0 and 1
  is_between_0_and_1 <- all(vec > 0 & vec < 1)

  # see if elements are strictly decreasing
  is_strictly_decreasing <- all(diff(vec) < 0)

  is_valid <- is_between_0_and_1 && is_strictly_decreasing

  return(is_valid)
}
