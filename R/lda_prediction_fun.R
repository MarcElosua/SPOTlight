#' This function processes the data and trains the LDA model
#'
#' @param lda_mod Object of class LDA_Gibbs.
#' @param spot_counts Object of class sparse matrix. The count matrix of the spots we want to predict with the classic GENESxSPOT
#' @param ncores Object of class integer, number of cores to use when parallelizing, by default it is 60% of the detected cores.
#' @param parallelize Object of class logical on wether to parallelize or not, by default TRUE.
#' @return matrix with the topic probability distribution for each spot.
#' @export
#' @examples
#'

lda_prediction <- function(lda_mod, spot_counts, ncores, parallelize=T){

  # Check variables
  if(is(se_obj)!="Seurat") {stop("ERROR: se_obj must be a Seurat object!")}
  # if(is(spot_counts)!="LDA_Gibbs") {stop("ERROR: lda_mod must be a LDA_Gibbs object!")}
  if(!is.integer(ncores)) {stop("ERROR: ncores must be an integer!")}
  if(!is.logical(parallelize)) {stop("ERROR: parallelize must be logical!")}

  #load required packages
  suppressMessages(require(topicmodels))
  suppressMessages(require(Matrix))
  suppressMessages(require(doSNOW))
  if(parallelize) suppressMessages(require(foreach))
  if(parallelize) suppressMessages(require(doParallel))

  spot_counts <- Matrix(sce_9cells_qc@assays$data$counts,sparse = T)
  # Only use genes detected by the model
  spot_counts <- spot_counts[lda_mod@terms,]

  lda_mod@terms %in%
  if(parallelize){
    # Detect number of cores and use 60% of them
    ncores <- round(parallel::detectCores() * 0.60)
    # Set up the backend
    cl <- parallel::makeCluster(ncores)
    # Register the backend
    doSNOW::registerDoSNOW(cl)
  }

  pred_start <- Sys.time()
  Sys.time()

  print('Running predictions')

  ## Set progress bar ##
  iterations <- nrow(spot_counts)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  prediction <- foreach(index=seq(1,nrow(spot_counts),10),
                        .combine = 'rbind',
                        .packages=c('topicmodels','Matrix'),
                        .options.snow = opts) %dopar% {

    test_spots_pred <- topicmodels::posterior(object = lda_mod,
                                              newdata = t(spot_counts[,index:(index+9)]))
    return(test_spots_pred$topics)

  }

  if(parallelize) parallel::stopCluster(cl)
  print(sprintf('Time to predict: %s', difftime(Sys.time(), pred_start, units = 'mins')))

  return(prediction)
}
