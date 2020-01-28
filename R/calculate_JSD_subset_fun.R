#' This function takes in 2 matrices of equal columns and calculates all pairwise JSD between them
#'
#' @param prediction: Object of class matrix with n ((n spots)) rows and k (n topics) columns where each column is a probabilistic array
#' @param syn_spots_profiles: Object of class matrix with m (n synthetic spots) rows and k (n topics) columns where each column is a probabilistic array
#' @param JSD_indices: Object of class list of m (n spots) elements with with each element having row I indices of syn_spots_profiles
#' @return this function returns a matrix with n rows and I columns the JSD values for the comparisons performed for each prediction
#' @export
#' @examples
#'
calculate_JSD_subset <- function(prediction, syn_spots_profiles, JSD_indices){

  # Check variables
  if( !( is.matrix(prediction) | is.data.frame(prediction)) ) {stop("ERROR: prediction must be a matrix object!")}
  if( !( is.matrix(syn_spots_profiles) | is.data.frame(syn_spots_profiles)) ) {stop("ERROR: syn_spots_profiles must be a matrix object!")}
  if( ! is.list(JSD_indices) ) {stop("ERROR: JSD_indices must be a list object!")}

  #### Initialize JS matrix ####
  mtrx_JSD_full <- matrix(nrow = nrow(prediction), ncol = max(lengths(JSD_indices)))
  print('Calculating Jensen-Shannon Divergence')
  pb_JSD <- txtProgressBar(min = 0, max = nrow(prediction), style = 3)

  ##### Calculate Jensen-Shannon divergence of the subset of the data
  for (i in 1:nrow(prediction)) {
    # print(i)
    for (ii in 1:length(JSD_indices[[i]])) {
      # print(sprintf('nested:%s',ii))
      x <- rbind(prediction[i,], syn_spots_profiles[JSD_indices[[i]][ii],])
      mtrx_JSD_full[i,ii] <- JSD(x, unit = "log2")
    }
    # update progress bar
    setTxtProgressBar(pb_JSD,i)
  }
  close(pb_JSD)
  return(mtrx_JSD_full)
}
