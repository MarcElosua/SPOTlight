#' This function processes the data and trains the LDA model
#'
#' @param prediction Object of class matrix with the predicted topic probability distributions for each spot. Output from the lda_prediction function.
#' @param syn_spots_ls list obtained from the function syn_spot_comb_topic_fun.R. The 1st element is a matrix with topic profiles of all the synthetic spots generated, the 2nd element is the composition of each synthetic spot.
#' @param top_dist Object of class integer, on how many top euclidean distance are we going to calculate the JSD.
#' @param top_JSD Object of class integer, how many of the top spots according JSD distance are we going to use to determine the composition.
#' @return matrix with the exact predicted composition of each spot.
#' @export
#' @examples
#'


syn_spot_assignment <- function(prediction, syn_spots_ls, top_dist=1000, top_JSD=15){

  # Check variables
  if( !( is.matrix(prediction) | is.data.frame(prediction)) ) {stop("ERROR: prediction must be a matrix object!")}
  if( !is.list(syn_spots_ls) ) {stop("ERROR: syn_spots_ls must be the list obtained from the function syn_spot_comb_topic_fun().")}
  if( !is.numeric(top_dist) ) {stop("ERROR: top_dist must be an integer!")}
  if( !is.numeric(top_JSD) ) {stop("ERROR: top_JSD must be an integer!")}


  #load required packages
  suppressMessages(require(Matrix))
  suppressMessages(require(pdist))
  suppressMessages(require(philentropy))

  syn_spots_profiles <- syn_spots_ls[[1]]
  syn_spots_metadata <- as.matrix(syn_spots_ls[[2]])
  syn_spots_metadata[is.na(syn_spots_metadata)] <- 0

  if (top_dist > nrow(syn_spots_profiles)) {
    warning(sprintf('top_dist cannot be larger than the total number of synthetic generated spots. Setting top_dist to %s', nrow(syn_spots_profiles)), sep = '\n')
    top_JSD <- top_dist
  }

  if (top_JSD > top_dist) {
    warning(sprintf('top_JSD cannot be larger than top_dist. Setting top_JSD to %s', top_dist), sep = '\n')
    top_JSD <- top_dist
  }

  ##### Calculate all pairwise euclidean distances between the predicted and simulated topic profiles #####
  dist <- pdist(X=prediction,Y=syn_spots_profiles)
  dist_mtrx <- as.matrix(dist)

  JSD_start <- Sys.time()

  ##### Get list with indices of best euclidean distance for each predictions #####
  JSD_indices <- top_n_predictions(dist_mtrx = dist_mtrx, n = top_dist)

  #### Calculate JSD for the subset of best predictions according to Euclidean distance #####
  mtrx_JSD_full <- suppressMessages(calculate_JSD_subset(prediction = prediction, syn_spots_profiles = syn_spots_profiles, JSD_indices = JSD_indices))

  quants_JSD <- round(quantile(matrixStats::rowMins(mtrx_JSD_full,na.rm = TRUE),c(0.25,0.5,0.75)),5)
  cat(sprintf("Quantiles of the JSD between the best synthetic spot profile and each spot's topic profile are - %s[%s-%s]", quants_JSD[[2]], quants_JSD[[1]], quants_JSD[[3]]), sep = '\n')

  ##### Get the index for each list from JSD_indices with the lowest JSD #####
  min15_error <- Rfast::rownth(x = mtrx_JSD_full, elems = rep(top_JSD, nrow(mtrx_JSD_full)), na.rm = TRUE)
  min_indices_JSD <- lapply(1:length(min15_error), function(i) which(mtrx_JSD_full[i,] <= min15_error[i]) )

  ##### Get Spot composition #####
  spot_composition_mtrx <- matrix(nrow = length(min_indices_JSD), ncol = ncol(syn_spots_metadata))
  colnames(spot_composition_mtrx) <- colnames(syn_spots_metadata)

  for (i in 1:nrow(spot_composition_mtrx)) {
    # Determine how many predictions we are adding since if there is only 1 we cannot do colmeans and we just need to assign it.
    best_comp <- syn_spots_metadata[JSD_indices[[i]][min_indices_JSD[[i]]],]

    if(is.null(nrow(best_comp))) { spot_composition_mtrx[i,] <- best_comp }
    else { spot_composition_mtrx[i,] <- round(colMeans(best_comp,na.rm = TRUE),0) }

  }; rm(i)

  return(spot_composition_mtrx)
}
