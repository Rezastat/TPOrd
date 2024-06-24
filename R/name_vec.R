name_vec <- function(vec, prenam, num_categories, model_type, type = c("est", "var"), ordered = F){

  vnam.orig <- names(vec)

  if (model_type == "Proportional_Odds") {
    mt_name = "po"
    vnam.new = paste0(mt_name, ".", type, ".", prenam, ".", vnam.orig)
    names(vec) <- vnam.new

  } else if (model_type == "Adjacent_category") {
    mt_name = "ac"
    vnam.new = paste0(mt_name, ".", type, ".", prenam, ".", vnam.orig)
    names(vec) <- vnam.new

  } else if (model_type == "Stopping_Ratio") {
    mt_name = "sr"
    vnam.new = paste0(mt_name, ".", type, ".", prenam, ".", vnam.orig)
    names(vec) <- vnam.new

  } else if (model_type == "Stereotype_Regression") {
    mt_name = "sg"
    vnam.new = paste0(mt_name, ".", type, ".", prenam, ".", vnam.orig)
    if (ordered == F) {
      vec <- c(vec[(num_categories-2):1], vec[(2*num_categories-3):(num_categories-1)], vec[(2*num_categories-2):length(vec)])
    }

    names(vec) <- vnam.new

  }
  return(vec)
}
