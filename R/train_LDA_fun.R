#' This function processes the data and trains the LDA model
#'
#' @param se_obj Object of class Seurat with the data of interest
#' @param clust_vr Name of the variable containing the cell clustering
#' @param cluster_markers_all Object of class dataframe obtained from the function Seurat::FindAllMarkers()
#' @param alpha Object of class "numeric"; initial value for alpha, by default 0.01.
#' @param verbose Object  of  class "integer".   If  a  positive  integer,  then  the  progress  is  reported  every "integer" verbose iterations. If 0 (default), no output is generated during model fitting
#' @param estimate.beta Object of class "logical"; controls if beta, the term distribution of the topics, isfixed, by default equals TRUE.
#' @param save Object of class "integer".  If a positive integer the estimated model is saved all verbose iterations. If 0 (default), no output is generated during model fitting.
#' @param keep Object of class "integer"; if a positive integer, the log-likelihood is saved every keep iterations, by default 100.
#' @param nstart Object of class "integer". Number of repeated random starts.
#' @param best Object of class "logical"; if TRUE only the model with the maximum (posterior) likelihood is returned, by default equals TRUE.
#' @param delta Object of class "numeric"; initial value for delta, by default equals 0.1.
#' @param iter Object of class "integer"; number of Gibbs iterations, by default equals 5000.
#' @param burnin Object of class "integer"; number of omitted Gibbs iterations at beginning, by default equals 0.
#' @param ... Other arguments passed on to LDA mod
#' @param thin Object of class "integer"; number of omitted in-between Gibbs iterations, by default equals iter.
#' @return A Latent Dirichlet Allocation list of model/s depending on if return best is T or F + a vector of the clusters of the cells used to train the model.
#' @export
#' @examples
#'

train_lda <- function(se_obj, clust_vr, cluster_markers_all, al=0.01, verbose=1, estimate.beta=TRUE, save=0, keep=100, nstart=1, best=TRUE, delta=0.1, iter=2000, burnin=0, thin=NULL, ...) {

  # Check variables
  if (is(se_obj) != "Seurat") stop("ERROR: se_obj must be a Seurat object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  if (!is.data.frame(cluster_markers_all)) stop("ERROR: cluster_markers_all must be a data frame object returned from Seurat::FindAllMarkers()!")
  if (!is.numeric(al)) stop("ERROR: al must be of class numeric!")
  if (al <= 0) stop("ERROR: al greater than 0!")
  if (verbose < 0) stop("ERROR: verbose must be an integer greater or equal than 0!")
  if (!is.logical(estimate.beta)) stop("ERROR: estimate.beta must be of class logical!")
  if (!is.numeric(save)) stop("ERROR: save must be of class integer!")
  if (!is.numeric(keep)) stop("ERROR: keep must be of class integer!")
  if (!is.numeric(nstart)) stop("ERROR: nstart must be of class integer!")
  if (!is.logical(best)) stop("ERROR: best must be of class logical!")
  if (!is.numeric(delta)) stop("ERROR: delta must be of class numeric!")
  if (!is.numeric(iter)) stop("ERROR: iter must be of class integer!")
  if (!is.numeric(burnin)) stop("ERROR: burnin must be of class integer!")
  if (!(is.numeric(thin) | is.null(thin))) stop("ERROR: thin must be of class integer!")

  # load required packages
  suppressMessages(require(Seurat))
  suppressMessages(require(topicmodels))
  suppressMessages(require(dplyr))
  suppressMessages(require(purrr))
  suppressMessages(require(tibble))
  suppressMessages(require(Matrix))

  se_obj$seurat_clusters <- droplevels(factor(se_obj@meta.data[, clust_vr]))

  #### Setting common parameters ####
  k <- nlevels(droplevels(factor(se_obj$seurat_clusters)))

  #### Get dataset ready ####
  se_lda_ready <- prep_seobj_topic_fun(se_obj = se_obj)

  # Select all marker genes for each cluster AND compute their Z score
  ntop <- max(table(cluster_markers_all$cluster))
  cluster_markers <- suppressMessages(cut_markers2(markers = cluster_markers_all, ntop = ntop))

  # Select unique markers from each cluster, if there are common markers between clusters lda model gets confused and classifies very different clusters as belonging to the same topic just because the seeding induced it!
  cluster_markers_uniq <- lapply(unique(cluster_markers$cluster), function(clust) {
    ls1 <- cluster_markers[cluster_markers$cluster == clust, "gene"]
    ls2 <- cluster_markers[cluster_markers$cluster != clust, "gene"]
    ls1_unique <- ls1[! ls1 %in% ls2]

    return(cluster_markers[cluster_markers$cluster == clust & cluster_markers$gene %in% ls1_unique, ])
  }) %>%
    bind_rows()

  # Set seedwords from top markers. Here we are setting the weights for each topic, the words that are weighted positively are those belonging to the list of top markers for a cluster.
  # In the seedgenes matrix each row represents a topic and each column represents a gene.

  # To the LDA model we need to pass a matrix with k rows and ngene columns, where each cell has the weight of that gene for that topic. The weight we're assigning is the logFC

  # initialize matrix
  seedgenes <- matrix(nrow = k, ncol = ncol(se_lda_ready), data = 0)
  colnames(seedgenes) = colnames(se_lda_ready)

  # Add seeds to model, if a cluster-topic has 0 unique markers its row will be set to all 0
  for (i in seq_len(k)) {
    clust_row <- cluster_markers_uniq$cluster == as.character(unique(se_obj@meta.data[, clust_vr])[[i]])
    seedgenes[i, cluster_markers_uniq[clust_row, "gene"]] = cluster_markers_uniq[clust_row, "logFC_z"]
    }

  # Verify that weights have been added
  # table(seedgenes != 0)

  #### LDA model ####
  # Set parameters
  control_LDA_Gibbs <- list(alpha = al, estimate.beta = estimate.beta,
                            verbose = verbose, prefix = tempfile(),
                            keep = keep, nstart = nstart, best = best,
                            delta = delta, iter = iter, burnin = burnin,
                            thin = thin, save = save,
                            seed = round(runif(nstart, min = 1, max = 1000)))

  # Train model
  s_gibbs_start <- Sys.time()
  print(s_gibbs_start)
  lda_mod <- topicmodels::LDA(se_lda_ready, k = k,
                 method = "Gibbs", seedwords = seedgenes, # Seedwords are only available with Gibbs sampling
                 control = control_LDA_Gibbs)
  print(sprintf("LDA seeded took: %s minutes",
                round(difftime(Sys.time(), s_gibbs_start, units = "mins"), 2))) # Takes ~10min

  if (is(lda_mod, "Gibbs_list")) {
    return(list(lda_mod@fitted, se_obj@meta.data[, clust_vr]))
  }  else {
    return(list(lda_mod, se_obj@meta.data[, clust_vr]))
  }
}
