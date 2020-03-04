#' This function processes the data and trains the LDA model
#'
#' @param lda_mod Object of class LDA_Gibbs.
#' @param spot_counts Object of class sparse matrix. The count matrix of the spots we want to predict as SPOTxGENES
#' @param ncores Object of class integer, number of cores to use when parallelizing, by default it is 60\% of the detected cores.
#' @return matrix with the topic probability distribution for each spot.
#' @export
#' @examples
#'

lda_prediction <- function(lda_mod, spot_counts, ncores) {

  # Check variables
  if (is(lda_mod)[[1]] != "LDA_Gibbs") stop("ERROR: lda_mod must be a LDA_Gibbs object!")
  if (!is.numeric(ncores)) stop("ERROR: ncores must be an integer!")

  #load required packages
  suppressMessages(require(topicmodels))
  suppressMessages(require(Matrix))
  suppressMessages(require(doSNOW))
  if (!is.null(ncores)) suppressMessages(require(foreach))
  if (!is.null(ncores)) suppressMessages(require(doParallel))

  # Set up cluster to parallelize
  if (!is.null(ncores)) {
    # Set up the backend
    cl <- parallel::makeCluster(ncores)
    # Register the backend
    doSNOW::registerDoSNOW(cl)
  }

  pred_start <- Sys.time()
  Sys.time()

  print("Running predictions")

  ## Set progress bar ##
  iterations <- length(seq(1, nrow(spot_counts), 10))
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  prediction <- foreach::foreach(index = seq(1, nrow(spot_counts), 10),
                        .combine = "rbind",
                        .options.snow = opts,
                        .packages = c("topicmodels", "Matrix", "dplyr")) %dopar% {

    index_end <- dplyr::if_else((index + 9) <= nrow(spot_counts),
                                as.double(index + 9), as.double(nrow(spot_counts)))

    test_spots_pred <- topicmodels::posterior(object = lda_mod,
                                              newdata = spot_counts[index:index_end, ])
    return(test_spots_pred$topics)

  }

  if (!is.null(ncores)) parallel::stopCluster(cl)
  cat(sprintf("\nTime to predict: %s minutes",
                round(difftime(Sys.time(), pred_start, units = "mins")), 2), sep = "\n")

  return(prediction)
}
