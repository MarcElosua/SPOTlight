#' This function processes the data and trains the LDA model
#'
#' @param prediction Object of class matrix with the predicted topic probability distributions for each spot. Output from the lda_prediction function.
#' @param syn_spots_ls list obtained from the function syn_spot_comb_topic_fun.R. The 1st element is a matrix with topic profiles of all the synthetic spots generated, the 2nd element is the composition of each synthetic spot.
#' @param top_dist Object of class integer, on how many top euclidean distance are we going to calculate the JSD.
#' @param top_jsd Object of class integer, how many of the top spots according JSD distance are we going to use to determine the composition.
#' @return matrix with the exact predicted composition of each spot.
#' @export
#' @examples
#'

syn_spot_assignment <- function(prediction, syn_spots_ls, top_dist=100, top_jsd=15) {

  # Check variables
  if (!(is.matrix(prediction) | is.data.frame(prediction))) stop("ERROR: prediction must be a matrix object!")
  if (!is.list(syn_spots_ls)) stop("ERROR: syn_spots_ls must be the list obtained from the function syn_spot_comb_topic_fun().")
  if (!is.numeric(top_dist)) stop("ERROR: top_dist must be an integer!")
  if (!is.numeric(top_jsd)) stop("ERROR: top_jsd must be an integer!")


  # load required packages
  suppressMessages(require(Matrix))
  suppressMessages(require(matrixStats))
  suppressMessages(require(pdist))
  suppressMessages(require(philentropy))

  # Extracting data
  syn_spots_profiles <- syn_spots_ls[[1]]
  syn_spots_metadata <- as.matrix(syn_spots_ls[[2]])
  syn_spots_metadata[is.na(syn_spots_metadata)] <- 0

  min_jsd_vec <- vector()

  ##### Get Spot composition #####
  spot_composition_mtrx <- matrix(nrow = nrow(prediction), ncol = ncol(syn_spots_metadata))
  colnames(spot_composition_mtrx) <- colnames(syn_spots_metadata)

  step <- 10
  iterations <- length(seq(1, nrow(prediction), step+1))
  pb <- txtProgressBar(min = 0, max = nrow(prediction), style = 3)

  # Iterate over the predictions 10 by 10 so as to not overload memory when doing pdist
  for (i in seq(1, nrow(prediction), step)) {

    #### Set index to nrow(prediction) in the last iteration ####
    index_end <- dplyr::if_else((i + step - 1) <= nrow(prediction),
                                as.double(i + step - 1),
                                as.double(nrow(prediction)))

    #### Subset the predictions to deconvolute ####
    prediction_subs <- prediction[i:index_end, ]

    ##### Calculate all pairwise euclidean distances between the predicted and simulated topic profiles #####
    dist <- pdist::pdist(X = prediction_subs,
                         Y = syn_spots_profiles)
    dist_mtrx <- as.matrix(dist)

    ##### Get list with indices of best euclidean distance for each predictions #####
    jsd_indices <- top_n_predictions(dist_mtrx = dist_mtrx,
                                     n = top_dist)

    #### Calculate JSD for the subset of best predictions according to Euclidean distance #####
    mtrx_jsd_full <- suppressMessages(calculate_jsd_subset(prediction = prediction_subs,
                                                           syn_spots_profiles = syn_spots_profiles,
                                                           jsd_indices = jsd_indices))

    #### Save the min jsd to get quantiles at the end ####
    min_jsd_vec <- c(min_jsd_vec, matrixStats::rowMins(mtrx_jsd_full, na.rm = TRUE))

    ##### Get the index for each list from JSD_indices with the lowest JSD #####
    min_jsd_error <- Rfast::rownth(x = mtrx_jsd_full,
                                   elems = rep(top_jsd, nrow(mtrx_jsd_full)),
                                   na.rm = TRUE)
    min_indices_jsd <- lapply(seq_len(length(min_jsd_error)), function(i) which(mtrx_jsd_full[i, ] <= min_jsd_error[i]))

    for (ii in seq_len(nrow(prediction_subs))) {
      # Determine how many predictions we are adding since if there is only 1 (top_jsd = 1) we cannot do colmeans when subsetting one row (it becomes a vector) and we just need to assign it.
      best_comp <- syn_spots_metadata[jsd_indices[[ii]][min_indices_jsd[[ii]]], ]

      if (is.null(nrow(best_comp))) {
        spot_composition_mtrx[i + ii - 1, ] <- best_comp
      } else {
        spot_composition_mtrx[i + ii - 1, ] <- round(colMeans(best_comp, na.rm = TRUE), 0)
      }

    }; rm(ii)
    # Update progress bar
    setTxtProgressBar(pb, i)
  }; rm(i)
  # Close progress bar
  close(pb)

  quants_jsd <- round(quantile(min_jsd_vec, c(0.25, 0.5, 0.75)), 5)
  print(sprintf("Quantiles of the JSD between the best synthetic spot profile and each spot's topic profile are - %s[%s-%s]",
                quants_jsd[[2]], quants_jsd[[1]], quants_jsd[[3]]))

  return(spot_composition_mtrx)
}
