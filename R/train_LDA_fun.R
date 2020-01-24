#' This function processes the data and trains the LDA model
#'
#' @param se_obj Object of class Seurat with the data of interest
#' @param clust_vr Name of the variable containing the cell clustering
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
#' @param thin Object of class "integer"; number of omitted in-between Gibbs iterations, by default equals iter.
#' @return A Latent Dirichlet Allocation list of model/s depending on if return best is T or F.
#' @export
#' @examples
#'

train_lda <- function(se_obj, clust_vr, al=0.01, estimate.beta=TRUE, save=0, keep=100, nstart=1, best=FALSE, delta=0.1, iter=5000, burnin=0, thin=NULL){

  # Check variables
  if(is(se_obj)!="Seurat") {stop("ERROR: se_obj must be a Seurat object!")}
  if(!is.character(clust_vr)){stop("ERROR: clust_vr must be a character string!")}
  if(!is.numeric(al)){stop("ERROR: al must be of class numeric!")}
  if(al <= 0){stop("ERROR: al greater than 0!")}
  if(!is.logical(estimate.beta)){stop("ERROR: estimate.beta must be of class logical!")}
  if(!is.integer(save)){stop("ERROR: save must be of class integer!")}
  if(!is.integer(keep)){stop("ERROR: keep must be of class integer!")}
  if(!is.integer(nstart)){stop("ERROR: nstart must be of class integer!")}
  if(!is.logical(best)){stop("ERROR: best must be of class logical!")}
  if(!is.numeric(delta)){stop("ERROR: delta must be of class numeric!")}
  if(!is.integer(iter)){stop("ERROR: iter must be of class integer!")}
  if(!is.integer(burnin)){stop("ERROR: burnin must be of class integer!")}
  if(!(is.integer(thin) | is.null(thin))){stop("ERROR: thin must be of class integer!")}

  #load required packages
  cat("Loading packages...", sep="\n")
  suppressMessages(require(Seurat))
  suppressMessages(require(topicmodels))
  suppressMessages(require(dplyr))
  suppressMessages(require(purrr))
  suppressMessages(require(Matrix))

  se_obj$seurat_clusters <- droplevels(se_obj@meta.data[,clust_vr])

  if(length(VariableFeatures(se_obj)) == 0) se_obj <- Seurat::FindVariableFeatures(object = se_obj, nfeatures = 5000)

  #### Extract the top marker genes from each cluster ####
  Idents(object = se_obj) <- se_obj$seurat_clusters
  cluster_markers_all <- Seurat::FindAllMarkers(object = se_obj, verbose = TRUE, only.pos = T)
  # saveRDS(object = cluster_markers_all,file = sprintf('%s/%s/cluster_markers_all_%s.RDS', an_01, robj_dir, ver))
  # cluster_markers_all <- readRDS(file = sprintf('%s/%s/cluster_markers_all_%s.RDS', an_01, robj_dir, ver))

  #### Combine marker genes and highest variable genes and subset genes ####
  keep_genes <- unique(c(VariableFeatures(se_obj), cluster_markers_all$gene))

  #### Get cell IDs to subset by cluster ####

  keep_ids <- lapply(split(se_obj@meta.data, se_obj@meta.data$seurat_clusters), function(subdf){
    # Determine n_sample, if the size of the group is < 100 use all the group, if not just use 100
    n_sample <- if_else(nrow(subdf) < 100, nrow(subdf), 100L)
    # Subset a random selection of that group and get the identifiers
    tmp_ds <- subdf[sample(1:nrow(subdf), n_sample),] %>%
      rownames_to_column('ID') %>%
      dplyr::pull(ID)
    return(tmp_ds)
  }) %>%
    purrr::flatten_chr() # flatten the list into a vector

  # keep_ids <- se_obj@meta.data %>%
  #   rownames_to_column('ID') %>%
  #   group_by(seurat_clusters) %>%
  #   sample_n(250) %>%
  #   ungroup() %>%
  #   dplyr::pull(ID)


  #### Subset seurat object ####
  se_obj <- se_obj[keep_genes,keep_ids]
  # saveRDS(se_obj, file = sprintf('%s/%s/seurat_umap_%s.RDS', an_00, robj_dir, id_comp))

  #### Setting common parameters ####
  k <- nlevels(droplevels(se_obj$seurat_clusters))
  nfeatures <- nrow(se_obj)

  #### Get dataset ready ####
  se_lda_ready <- prep_seobj_topic_fun(se_obj = se_obj)
  # se_obj <- FindVariableFeatures(object = se_obj, nfeatures = 5000)
  # se_obj <- SCTransform(se_obj, verbose = T,variable.features.n = 5000)

  # Select top 100 genes for each cluster
  cluster_markers <- cut_markers2(clusters = levels(cluster_markers_all$cluster), markers = cluster_markers_all, ntop = 100)

  # Select unique markers from each cluster, if there are common markers between clusters lda model gets confused and classifies very different clusters as belonging to the same topic just because the seeding induced it!
  cluster_markers_uniq <- lapply(unique(cluster_markers$cluster), function(clust){
    ls1 <- cluster_markers[cluster_markers$cluster == clust,'gene']
    ls2 <- cluster_markers[cluster_markers$cluster != clust,'gene']
    ls1_unique <- ls1[!ls1 %in% ls2]

    return(cluster_markers[cluster_markers$cluster == clust & cluster_markers$gene %in% ls1_unique,])
  }) %>%
    bind_rows()

  # Set seedwords from top markers. Here we are setting the weights for each topic, the words that are weighted positively are those belonging to the list of top markers for a cluster.
  # In the seedgenes matrix each row represents a topic and each column represents a gene.

  # To the LDA model we need to pass a matrix with k rows and ngene columns, where each cell has the weight of that gene for that topic. The weight we're assigning is the logFC

  # initialize matrix
  seedgenes <- matrix(nrow=k, ncol=ncol(se_lda_ready), data=0)
  colnames(seedgenes) = colnames(se_lda_ready)


  for (i in 1:k) { seedgenes[i,cluster_markers_uniq[cluster_markers_uniq$cluster == cluster_markers_uniq$cluster[[i]],'gene']] = cluster_markers_uniq[cluster_markers_uniq$cluster == cluster_markers_uniq$cluster[[i]],'logFC_z']; print(i) }

  # Verify that weights have been added
  table(seedgenes != 0)

  #### LDA model ####
  # Set parameters
  control_LDA_Gibbs <- list(alpha = al, estimate.beta = TRUE,
                            verbose = 0, prefix = tempfile(), save = 0, keep = 0,
                            seed = as.integer(Sys.time()), nstart = 1, best = F,
                            delta = 0.1, iter = 5000, burnin = 3000, thin = 100)

  # Train model
  s_gibbs_seed <- Sys.time()
  print(s_gibbs_seed)
  set.seed(123)
  lda_mod <- LDA(se_lda_ready, k=k,
                 method='Gibbs', seedwords=seedgenes, # Seedwords are only available with Gibbs sampling
                 control=control_LDA_Gibbs)
  print(sprintf('LDA seeded took: %s', difftime(Sys.time(), s_gibbs_seed, units = 'mins'))) # Takes ~10min

  if(is(lda_mod) == "Gibbs_list") return(lda_mod@fitted) else return(list(lda_mod))
}
