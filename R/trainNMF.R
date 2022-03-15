#' @name trainNMF
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
#' @param assay if the object is of Class \code{Seurat}, character string
#'   specifying the assay from which to extract the expression matrix.
#'     By default "RNA".
#' @param slot if the object is of Class \code{Seurat}, character string
#'   specifying the slot from which to extract the expression matrix. If the
#'   object is of class \code{SingleCellExperiment} indicates matrix to use.
#'   By default "counts".
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
#' # mock up some single-cell, mixture & marker data
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

#' @rdname trainNMF
#' @importFrom SingleCellExperiment colLabels
#' @export
setMethod("trainNMF",
    c("SingleCellExperiment", "ANY"),
    function(x, y, ...,
        assay = "counts",
        groups = colLabels(x, onAbsence = "error"))
    {
        # Check necessary packages are installed and if not STOP
        .test_installed("SummarizedExperiment")
        trainNMF(as.matrix(SummarizedExperiment::assay(x, assay)),
            y, groups, ...)
    })

#' @rdname trainNMF
#' @export
setMethod("trainNMF",
    c("ANY", "SingleCellExperiment"),
    function(x, y, ...,
        assay = "counts") {
        # Check necessary packages are installed and if not STOP
        .test_installed("SummarizedExperiment")
        trainNMF(x, as.matrix(SummarizedExperiment::assay(y, assay)), ...)
    })

#' @rdname trainNMF
#' @importFrom SeuratObject Idents GetAssayData
#' @export
setMethod("trainNMF",
    c("Seurat", "ANY"),
    function(x, y, ...,
        slot = "counts",
        assay = "RNA",
        groups = Idents(x)) {
        .test_installed("SeuratObject")
        trainNMF(GetAssayData(x, slot, assay), y, groups, ...)
    })

#' @rdname trainNMF
#' @importFrom SeuratObject GetAssayData
#' @export
setMethod("trainNMF",
    c("ANY", "Seurat"),
    function(x, y, ...,
        slot = "counts",
        assay = "RNA") {
        .test_installed("SeuratObject")
        trainNMF(x, GetAssayData(y, slot, assay), ...)
    })

#' @rdname trainNMF
#' @export
setMethod("trainNMF",
    c("ANY", "dgCMatrix"),
    function(x, y, ...,
        slot = "counts",
        assay = "RNA") {
        trainNMF(x, as.matrix(y), ...)
    })

#' @rdname trainNMF
#' @export
setMethod("trainNMF",
    c("dgCMatrix", "ANY"),
    function(x, y, ...,
        slot = "counts",
        assay = "RNA") {
        trainNMF(as.matrix(x), y, ...)
    })

#' @rdname trainNMF
#' @export
setMethod("trainNMF",
    c("ANY", "DelayedMatrix"),
    function(x, y, ...,
        slot = "counts",
        assay = "RNA") {
        trainNMF(x, as.matrix(y), ...)
    })

#' @rdname trainNMF
#' @export
setMethod("trainNMF",
    c("DelayedMatrix", "ANY"),
    function(x, y, ...,
        slot = "counts",
        assay = "RNA") {
        trainNMF(as.matrix(x), y, ...)
    })

#' @rdname trainNMF
#' @export
setMethod("trainNMF",
    c("ANY", "ANY"),
    function(x, y, ...) {
        stop("See ?trainNMF for valid x & y inputs")
    })

#' @rdname trainNMF
#' @importFrom Matrix rowSums
#' @importFrom NMF nmf nmfModel
#' @export
setMethod("trainNMF",
    c("matrix", "matrix"),
    function(x, y,
        groups,
        mgs,
        n_top = NULL,
        gene_id = "gene",
        group_id = "cluster",
        weight_id = "weight",
        hvg = NULL,
        model = c("ns", "std"),
        scale = TRUE,
        verbose = TRUE,
        ...) {
        # check validity of input arguments
        model <- match.arg(model)
    
    # select genes in mgs or hvg
    if (!is.null(hvg)) {
        # Select union of genes between markers and HVG
        mod_genes <- union(unique(mgs[, gene_id]), hvg)
    } else {
        # Select genes from the marker genes
        mod_genes <- unique(mgs[, gene_id])
    }
    
    # Select intersection between interest and present in x (sce) & y (spe)
    mod_genes <- intersect(mod_genes, intersect(rownames(x), rownames(y)))
    
    # drop features that are undetected
    # in single-cell and/or mixture data
    x <- .filter(x[mod_genes, ], y[mod_genes, ])
    mgs <- mgs[mgs[[gene_id]] %in% rownames(x), ]
    
    # scale to unit variance (optional)
    if (scale) {
        if (verbose) message("Scaling count matrix")
        x <- .scale_uv(x)
    }
    
    # capture start time
    t0 <- Sys.time()
    
    # set model rank to number of groups
    rank <- length(unique(groups))
    
    # get seeding matrices (optional)
    seed <- if (TRUE) {
        if (verbose) message("Seeding initial matrices")
        hw <- .init_nmf(x, groups, mgs, n_top, gene_id, group_id, weight_id)
        nmfModel(W = hw$W, H = hw$H, model = paste0("NMF", model))
    }
    
    # train NMF model
    if (verbose) message("Training NMF model")
    mod <- nmf(x, rank, paste0(model, "NMF"), seed, ...)
    
    # capture stop time
    t1 <- Sys.time()
    
    # print runtimes
    if (verbose) {
        dt <- round(difftime(t1, t0, units = "mins"), 2)
        message("Time for training: ", dt, "min")
    }
    
    # get topic profiles
    topic <- .topic_profiles(mod, groups)
    
    return(list("mod" = mod, "topic" = topic))
})