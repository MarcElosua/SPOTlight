#' This function gets a distance matrix and returns a list with the indices of the lowest (best) prediction indices
#'
#' @param dist_mtrx Object of class matrix. Distance matrix to evaluate
#' @param n number of best predictions to return
#' @return list of nrows(dist_mtrx) elements with ~n elements each
#' @export
#' @examples
#'

top_n_predictions <- function(dist_mtrx, n) {

  if (!is.matrix(dist_mtrx)) stop("ERROR: dist_mtrx must be a matrix object!")

  # check that n isn't bigger than nrow(dist_mtrx), that will throw it off
  n <- dplyr::if_else(n < ncol(dist_mtrx), n, as.double(ncol(dist_mtrx)))
  min_error <- Rfast::rownth(x = dist_mtrx, elems = rep(n, nrow(dist_mtrx)))

  # Get indices over which to calculate JD
  jd_indices <- lapply(seq_len(length(min_error)), function(i) {
    which(dist_mtrx[i, ] <= min_error[i])
  })
  return(jd_indices)
}
