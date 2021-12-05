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
#' @param hvg character vector containing hvg to include in the model. 
#'   By default NULL.
#' @param gene_id,group_id,weight_id character specifying the column in \code{mgs}
#'   containing gene identifiers, group labels and weights, respectively.
#' @param scale logical specifying whether to scale single-cell counts to unit
#'   variance. This gives the user the option to normalize the data beforehand
#'   as you see fit (CPM, FPKM, ...) when passing a matrix or specifying the
#'   slot from where to extract the count data.
#' @param min_prop scalar in [0,1] setting the minimum contribution
#'   expected from a cell type in \code{x} to observations in \code{y}.
#'   By default 0.
#' @param verbose logical. Should information on progress be reported?
#' 
#' @return a numeric matrix with rows corresponding to samples and columns to groups
#' 
#' @author Marc Elosua-Bayes & Helena L. Crowell
#' 
#' @details SPOTlight uses a Non-Negative Matrix Factorization approach to learn
#'   which genes are important for each cell type. In order to drive the
#'   factorization and give more importance to cell type marker genes we
#'   previously compute them and use them to initialize the basis matrix. This
#'   initialized matrices will then be used to carry out the factorization with 
#'   the single cell expression data. Once the model has learn the topic
#'   profiles for each cell type we use non-negative least squares (NNLS) to
#'   obtain the topic contributions to each spot. Lastly, NNLS is again used to 
#'   obtain the proportion of each cell type for each spot by finding the
#'   fitting the single-cell topic profiles to the spots topic contributions.
#' 
# TODO grab ideas from vignette/unit tests
#' @examples 
#' # Get sc Brain data from 
#' # TENxBrainData
#' # http://127.0.0.1:25293/library/TENxBrainData/doc/TENxBrainData.html
#' # ExperimentHub
#' # https://bioconductor.org/packages/release/bioc/vignettes/ExperimentHub/inst/doc/ExperimentHub.html#using-experimenthub-to-retrieve-data
#' library(TENxBrainData)
#' library(ExperimentHub)
#' sc <- TENxBrainData20k()
#' 
#' # Get Visium data from TENxVisiumData
#' # https://bioconductor.org/packages/release/data/experiment/vignettes/TENxVisiumData/inst/doc/vignette.html
#' library(TENxVisiumData)
#'TODO 
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
#' @importFrom SeuratObject Idents GetAssayData
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
#' @importFrom SeuratObject GetAssayData
#' @export
setMethod("SPOTlight",
    c("ANY", "Seurat"),
    function(x, y, ..., 
        slot = "counts", 
        assay = "RNA") 
    {
        SPOTlight(x, GetAssayData(y, slot, assay), ...)
    })

#' @rdname SPOTlight
#' @export
setMethod("SPOTlight",
  c("ANY", "dgCMatrix"),
  function(x, y, ..., 
    slot = "counts", 
    assay = "RNA") 
  {
    SPOTlight(x, as.matrix(y), ...)
  })

#' @rdname SPOTlight
#' @export
setMethod("SPOTlight",
  c("dgCMatrix", "ANY"),
  function(x, y, ..., 
    slot = "counts", 
    assay = "RNA") 
  {
    SPOTlight(as.matrix(x), y, ...)
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
        hvg = NULL,
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
        
        # train NMF model
        mod <- .train_nmf(x, y, groups, mgs, n_top, gene_id, group_id, weight_id, hvg, model, scale, verbose)
        
        # get topic profiles
        ref <- .topic_profiles(mod, groups)
        
        # perform deconvolution
        res <- .deconvolute(y, mod, ref, scale, min_prop, verbose)
        
        # return list of NMF model & deconvolution matrix
        list(res, mod)
    })
