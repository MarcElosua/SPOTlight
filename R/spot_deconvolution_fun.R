#' This function processes the data and trains the LDA model
#'
#' @param lda_mod Object of class LDA_Gibbs.
#' @param train_cell_clust Object of class vector with cluster of the cells used to train the model.
#' @param clust_vr Object of class character. Name of the variable containing the cell clustering.
#' @param spot_counts Object of class sparse matrix. Count matrix of the spots to deconvolute.
#' @param verbose Object of class Logical determining if progress should be reported or not (TRUE by default).
#' @param ncores Object of class integer, How many cores to use to parallelize the process, if NULL i the process will be run sequentially.
#' @param top_dist Object of class integer, on how many top euclidean distance are we going to calculate the JSD.
#' @param top_jsd Object of class integer, how many of the top spots according JSD distance are we going to use to determine the composition.
#' @return matrix with the exact predicted composition of each spot.
#' @export
#' @examples
#'

spot_deconvolution <- function(lda_mod, train_cell_clust, clust_vr, spot_counts,
                               verbose = TRUE, ncores = NULL,
                               top_dist = 1000, top_jsd = 10) {

  # Check variables
  if (!is(lda_mod, "LDA_Gibbs")) stop("ERROR: lda_mod must be an LDA_Gibbs object!")
  if (! is(train_cell_clust, "vector")) stop("ERROR: se_obj must be a vector/list object!")
  if (!is.character(clust_vr))stop("ERROR: clust_vr must be a character string!")
  if (!is.logical(verbose))stop("ERROR: verbose must be a logical object!")
  if (!(is.numeric(ncores) | is.null(ncores))) stop("ERROR: ncores must be an integer!")
  if (!is.numeric(top_dist)) stop("ERROR: top_dist must be an integer!")
  if (!is.numeric(top_jsd)) stop("ERROR: top_jsd must be an integer!")


  #load required packages
  suppressMessages(require(pdist))
  suppressMessages(require(philentropy))
  suppressMessages(require(Matrix))
  suppressMessages(require(arrangements))
  suppressMessages(require(progress))
  suppressMessages(require(purrr))

  # Create all synthetic spot combinations
  syn_spots_ls <- syn_spot_comb_topic(lda_mod = lda_mod,
                                      train_cell_clust = train_cell_clust,
                                      clust_vr = clust_vr,
                                      verbose = verbose)

  # Predict topic profiles of spatial spots
  if (verbose) print("Predict topic profiles of spatial spots")
  prediction <- lda_prediction(lda_mod = lda_mod,
                               spot_counts = spot_counts,
                               ncores = ncores)

  # Perform deconvolution of the spatial spots
  if (verbose) print("Perform deconvolution of the spatial spots")
  spot_deconv <- syn_spot_assignment(prediction = prediction,
                                     syn_spots_ls = syn_spots_ls,
                                     top_dist = top_dist,
                                     top_jsd = top_jsd)

  return(spot_deconv)

}
