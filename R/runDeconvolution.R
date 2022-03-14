#' @rdname trainNMF
#' @title train NMF model
#'
#' @aliases trainNMF
#'
#' @description This is the training function used by SPOTLight. This function
#'   takes in single cell expression data, trains the model and learns topic
#'    profiles for each cell type
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
#' @param gene_id,group_id,weight_id character specifying the column
#'   in \code{mgs} containing gene identifiers, group labels and weights,
#'   respectively.
#' @param scale logical specifying whether to scale single-cell counts to unit
#'   variance. This gives the user the option to normalize the data beforehand
#'   as you see fit (CPM, FPKM, ...) when passing a matrix or specifying the
#'   slot from where to extract the count data.
#' @param n_top integer scalar specifying the number of markers to select per
#'  group. By default NULL uses all the marker genes to initialize the model.
#' @param model character string indicating which model to use when running NMF.
#' Either "ns" (default) or "std".
#' @param verbose logical. Should information on progress be reported?
#' @param ... additional parameters.
#'
#'
#' @return base a list where the first element is an \code{NMFfit} object and
#'   the second is a matrix contatining the topic profiles learnt.
#'
#' @author Marc Elosua Bayes & Helena L Crowell
#'
#' @examples
#' set.seed(321)
#' mock up some single-cell, mixture & marker data
#' sce <- mockSC(ng = 200, nc = 10, nt = 3)
#' spe <- mockSP(sce)
#' mgs <- getMGS(sce)
#' 
#' res <- trainNMF(
#'     x = sce,
#'     y = spe,
#'     groups = sce$type,
#'     mgs = mgs,
#'     weight_id = "weight",
#'     group_id = "type",
#'     gene_id = "gene")
#' # Get NMF model
#' res[["mod"]]
#' # Get topic profiles
#' res[["topic"]]
NULL

#' @rdname runDeconvolution
#' @export
setMethod("runDeconvolution", "SingleCellExperiment",
    function(x, ..., assay = "counts") {
        # Check necessary packages are installed and if not STOP
        .test_installed("SummarizedExperiment")
        runDeconvolution(as.matrix(SummarizedExperiment::assay(x, assay)), ...)
    })

#' @rdname runDeconvolution
#' @export
setMethod("runDeconvolution", "Seurat",
    function(x, slot = "counts", assay = "RNA", ...) {
        runDeconvolution(GetAssayData(x, slot, assay), ...)
    })

#' @rdname runDeconvolution
#' @export
setMethod("runDeconvolution", "dgCMatrix",
    function(x, ...) {
        runDeconvolution(as.matrix(x), ...)
    })

#' @rdname runDeconvolution
#' @export
setMethod("runDeconvolution", "DelayedMatrix",
    function(x, ...) {
        runDeconvolution(as.matrix(x), ...)
    })

#' @rdname runDeconvolution
#' @export
setMethod("runDeconvolution", "ANY",
    function(x, ...) {
        stop("See ?runDeconvolution for valid x & y inputs")
    })

#' @importFrom nnls nnls
#' #' @rdname runDeconvolution
#' @export
setMethod("runDeconvolution", "matrix",
    function(x, mod,
        ref,
        scale = TRUE,
        min_prop = 0.01,
        verbose = TRUE) {
        
        # Get topic profiles for mixtures
        mat <- .pred_prop(x, mod, scale)
        
        if (verbose) message("Deconvoluting mixture data")
        
        res <- vapply(seq_len(ncol(mat)), function(i) {
            pred <- nnls::nnls(ref, mat[, i])
            prop <- prop.table(pred$x)
            # drop groups that fall below 'min_prop' & update
            prop[prop < min_prop] <- 0
            prop <- prop.table(prop)
            # compute residual sum of squares
            ss <- sum(mat[, i]^2)
            # compute percentage of unexplained residuals
            err <- pred$deviance / ss
            c(prop, err)
        }, numeric(ncol(ref) + 1))
        
        # set dimension names
        rownames(res) <- c(dimnames(mod)[[3]], "res_ss")
        colnames(res) <- colnames(mat)
        
        # Separate residuals from proportions
        # Extract residuals
        err <- res["res_ss", ]
        # Extract only deconvolution matrices
        res <- res[-nrow(res), ]
        
        return(list("mat" = t(res), "res_ss" = err))
})
