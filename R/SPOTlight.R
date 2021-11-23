#' @name SPOTlight
#' @title Deconvolution of mixture using single-cell data
#' 
#' @description This is the backbone function which takes in single cell
#'   expression data to deconvolute spatial transcriptomics spots.
#' 
#' @param x,y single-cell and mixture dataset, respectively. Can be a 
#'   numeric matrix, \code{SingleCellExperiment} or \code{SeuratObjecy}.
#' @param groups vector of group labels for cells in \code{x}. 
#'   When \code{x} is a \code{SingleCellExperiment} or \code{SeuratObject}, 
#'   defaults to \code{colLabels} and \code{Idents(x)}, respectively.
#' @param mgs \code{data.frame} or \code{DataFrame} of marker genes. 
#'   Must contain columns holding gene identifiers, group labels and
#'   the weight (e.g., logFC, -log(p-value) a feature has in a given group. 
#' @param gene_id,group_id,weight_id character specifying the column in \code{mgs}
#'   containing gene identifiers, group labels and weights, respectively.
#' @param scale logical specifying whether to
#'   scale single-cell counts to unit variance.
#' @param min_prop scalar in [0,1] setting the minimum contribution
#'   expected from a cell type in \code{x} to observations in \code{y}.
#' @param verbose logical. Should information on progress be reported?
#' 
#' @return a numeric matrix with rows corresponding to samples and columns to groups
#' 
#' @author Marc Elosua-Bayes & Helena L. Crowell
#' 
#' @examples 
#' print("todo")
NULL

#' @rdname SPOTlight
#' @importFrom SingleCellExperiment colLabels
#' @importFrom SummarizedExperiment assay
#' @export
setMethod("SPOTlight", 
    c("SingleCellExperiment", "ANY"),
    function(x, y, ..., 
        assay = "counts", 
        groups = colLabels(x, onAbsence = "error")) 
    {
        SPOTlight(as.matrix(assay(x, assay)), y, groups, ...)
    })

#' @rdname SPOTlight
#' @importFrom SummarizedExperiment assay
#' @export
setMethod("SPOTlight", 
    c("ANY", "SingleCellExperiment"),
    function(x, y, ..., 
        assay = "counts") 
    {
        SPOTlight(x, as.matrix(assay(y, assay)), ...)
    })

#' @rdname SPOTlight
#' @importFrom Seurat Idents GetAssayData
#' @export
setMethod("SPOTlight",
    c("Seurat", "ANY"),
    function(x, y, ..., 
        slot = "counts", 
        assay = "RNA", 
        groups = Idents(x)) 
    {
        SPOTlight(GetAssayData(x, slot, assay), y, groups, ...)
    })

#' @rdname SPOTlight
#' @importFrom Seurat Idents GetAssayData
#' @export
setMethod("SPOTlight",
    c("ANY", "Seurat"),
    function(x, y, ..., 
        slot = "counts", 
        assay = "RNA") 
    {
        SPOTlight(x, GetAssayData(y, slot, assay), ...)
    })

setMethod("SPOTlight", 
    c("ANY", "ANY"), 
    function(x, y, ...) 
    {
        stop("...")
    })

#' @rdname SPOTlight
#' @export
setMethod("SPOTlight", 
    c("matrix", "matrix"),
    function(x, y, 
        groups,
        # markers
        mgs,
        n_top = NULL,
        gene_id = "gene",
        group_id = "cluster",
        weight_id = "weight",
        # NMF
        scale = TRUE,
        model = c("ns", "std"),
        # deconvolution
        min_prop = 0.01,
        # other
        verbose = TRUE)
    {
        # check validity if input arguments
        model <- match.arg(model)
        if (is.null(n_top)) 
            n_top <- max(table(mgs[[group_id]]))
        ids <- c(gene_id, group_id, weight_id)
        stopifnot(
            is.numeric(x), is.numeric(y),
            is.character(ids), length(ids) == 3, ids %in% names(mgs),
            length(groups) == ncol(x), groups %in% mgs[[group_id]],
            is.numeric(min_prop), length(min_prop) == 1,
            min_prop >= 0, min_prop <= 1,
            is.logical(scale), length(scale) == 1,
            is.logical(verbose), length(verbose) == 1)

        # # downsample scRNA-seq to select gene set 
        # # and number of cells to train the model  
        # sub <- .downsample_sce(sce, cells, genes)

        # train NMF model
        mod <- .train_nmf(x, y, groups, mgs, n_top, gene_id, group_id, weight_id, model, scale, verbose)
        
        # get topic profiles
        ref <- .topic_profiles(mod, groups)
        
        # perform deconvolution
        res <- .deconvolute(y, mod, ref, scale, min_prop, verbose)
        
        # return list of NMF model & deconvolution matrix
        list(mod, res)
    })