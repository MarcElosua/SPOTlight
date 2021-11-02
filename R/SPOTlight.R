#' @rdname SPOTlight
#' @title Deconvolution of mixture using single-cell data
#' 
#' @Description ...
#' 
#' @param sce \code{\link{SingleCellExperiment} of single-cell data
#' @param spe \code{\link{SpatialExperiment} of mixture data
#' 
#' @examples 
#' print("todo")
#' 
#' @author Marc Elosua Bayes
#' 
#' @importFrom methods is
#' @importFrom SingleCellExperiment colLabels colLabels<-
#' 
#' @export

SPOTlight <- function(
    sce,
    spe,
    groups = colLabels(sce, onAbsence = "error"),
    # filtering
    n_cells = 100,
    n_genes = 3e3,
    # markers
    mgs,
    n_top = NULL,
    gene_id = "gene",
    group_id = "cluster",
    weight_id = "weight",
    # NMF
    scale = TRUE,
    seed_model = c("ns", "std"),
    NMF_model
    # deconvolution
    min_prop = 0.01,
    # general
    assay = "counts",
    verbose = TRUE)
{
    # check validity if input arguments
    stopifnot(
        is(sce, "SingleCellExperiment"),
        #is(spe, "SpatialExperiment"),
        is.logical(scale), length(scale) == 1,
        is.logical(verbose), length(verbose) == 1)
    
    if (length(groups) == 1) {
        stopifnot(
            is.character(groups), 
            !is.null(sce[[groups]]))
    } else if (length(groups) == ncol(sce)) {
        stopifnot(groups %in% mgs[[group_id]])
    } else {
        stop("...")
    }
    
    # assure 'colLabels' are set
    if (is.null(colLabels(sce)))
        colLabels(sce) <- groups
    
    # downsample scRNA-seq to select gene set 
    # and number of cells to train the model  
    sub <- .downsample_sce(sce, n_cells, n_genes)
    
    # train NMF model
    mod <- .train_nmf(sub, spe, mgs, scale = scale, verbose = verbose)
    
    # get topic profiles
    ref <- .topic_profiles(sub, mod)
    
    # perform deconvolution
    res <- .deconvolute(mod, spe, scale, ref, min_prop)
}
