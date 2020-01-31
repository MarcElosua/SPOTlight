#' This function takes in 2 matrices of equal columns and calculates all pairwise JSD between them
#'
#' @param prediction: Object of class matrix with n ((n spots)) rows and k (n topics) columns where each column is a probabilistic array
#' @param syn_spots_profiles: Object of class matrix with m (n synthetic spots) rows and k (n topics) columns where each column is a probabilistic array
#' @param jsd_indices: Object of class list of m (n spots) elements with with each element having row I indices of syn_spots_profiles
#' @return this function returns a matrix with n rows and I columns the JSD values for the comparisons performed for each prediction
#' @export
#' @examples
#'
calculate_jsd_subset <- function(prediction, syn_spots_profiles, jsd_indices) {

  # Check variables
  if (!(is.matrix(prediction) | is.data.frame(prediction))) stop("ERROR: prediction must be a matrix object!")
  if (!(is.matrix(syn_spots_profiles) | is.data.frame(syn_spots_profiles))) stop("ERROR: syn_spots_profiles must be a matrix object!")
  if (!is.list(jsd_indices)) stop("ERROR: jsd_indices must be a list object!")

  # load required packages
  suppressMessages(require(philentropy))

  #### Initialize JS matrix ####
  mtrx_jsd_full <- matrix(nrow = nrow(prediction),
                          ncol = max(lengths(jsd_indices)))
  print("Calculating Jensen-Shannon Divergence")
  pb_jsd <- txtProgressBar(min = 0, max = nrow(prediction), style = 3)

  ##### Calculate Jensen-Shannon divergence of the subset of the data
  for (i in seq_len(nrow(prediction))) {
    for (ii in seq_len(length(jsd_indices[[i]]))) {
      x <- rbind(prediction[i, ], syn_spots_profiles[jsd_indices[[i]][ii], ])
      mtrx_jsd_full[i, ii] <- suppressMessages(JSD(x, unit = "log2"))
    }
    # update progress bar
    setTxtProgressBar(pb_jsd, i)
  }
  close(pb_jsd)
  return(mtrx_jsd_full)
}
