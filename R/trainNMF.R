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
#' @param groups character vector of group labels for cells in \code{x}.
#'   When \code{x} is a \code{SingleCellExperiment} or \code{SeuratObject},
#'   defaults to \code{colLabels(x)} and \code{Idents(x)}, respectively.
#'   Make sure groups is not a Factor.
#' @param mgs \code{data.frame} or \code{DataFrame} of marker genes.
#'   Must contain columns holding gene identifiers, group labels and
#'   the weight (e.g., logFC, -log(p-value) a feature has in a given group.
#' @param pnmf character vector specifying which from which package to grab
#'   \code{nmf()}. It can be either \code{RcppML} (default) or \code{NMF}.
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
#' @param model character string indicating which model to use when running nmf
    #'   from the NMF. Either "ns" (default) or "std".
#' @param assay_sc,assay_sp if the object is of Class \code{Seurat}, character string
#'   specifying the assay from which to extract the expression matrix.
#'   By default "RNA" and "Spatial".
#' @param slot_sc,slot_sp if the object is of Class \code{Seurat}, character string
#'   specifying the slot from which to extract the expression matrix. If the
#'   object is of class \code{SingleCellExperiment} indicates matrix to use.
#'   By default "counts".
#' @param L1 LASSO penalties in the range (0, 1], single value or array of
#'   length two for c(w, h). See ?RcppML::nmf() for more info.
#' @param L2 Ridge penalties greater than zero, single value or array of length
#'   two for c(w, h). See ?RcppML::nmf() for more info.
#' @param tol tolerance of the fit ?RcppML::nmf() for more info.
#' @param verbose logical. Should information on progress be reported?
#' @param ... additional parameters.
#'
#'
#' @return a list where the first element is a list with the NMF model
#'   information (see ?RcppML::nmf return value) and the second is a matrix
#'    containing the topic profiles learnt per cell type.
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

#' @importFrom Matrix rowSums Matrix
#' @export
trainNMF <- function(
    x,
    y,
    groups = NULL,
    mgs,
    pnmf = c("RcppML", "NMF"),
    n_top = NULL,
    gene_id = "gene",
    group_id = "cluster",
    weight_id = "weight",
    hvg = NULL,
    model = c("ns", "std"),
    scale = TRUE,
    verbose = TRUE,
    L1_nmf = 0,
    L2_nmf = 0,
    tol = 1e-05,
    assay_sc = "RNA",
    slot_sc = "counts",
    assay_sp = "Spatial",
    slot_sp = "counts",
    ...) {
    # check validity of input arguments
    pnmf <- match.arg(pnmf)
    model <- match.arg(model)
    
    if (is.null(n_top))
        n_top <- max(table(mgs[[group_id]]))
    ids <- c(gene_id, group_id, weight_id)
    
    stopifnot(
        is.numeric(x) | is(x, "dgCMatrix") |
            is(x, "Seurat") | is(x, "SingleCellExperiment") |
            is(x, "DelayedMatrix"), 
        is.numeric(y) | is(y, "dgCMatrix") |
            is(y, "Seurat") | is(y, "SingleCellExperiment") |
            is(y, "DelayedMatrix"),
        is.character(ids), length(ids) == 3, ids %in% names(mgs),
        is.null(groups) | length(groups) == ncol(x),
        is.logical(scale), length(scale) == 1,
        is.logical(verbose), length(verbose) == 1,
        is.numeric(L1_nmf), length(L1_nmf) < 3,
        is.numeric(L2_nmf), length(L2_nmf) < 3,
        is.numeric(tol), length(tol) == 1)
    
    # Set groups if x is SCE or SE and groups is NULL 
    if (is.null(groups))
        groups <- .set_groups_if_null(x)
    
    groups <- as.character(groups)
    
    # Stop if at least one of the groups doesn't have marker genes
    stopifnot(groups %in% mgs[[group_id]])
    
    # Extract expression matrices for x and y
    # TODO check this step
    if (!is.matrix(x) & !is(x, "dgCMatrix"))
        x <- .extract_counts(x, assay_sc, slot_sc)
    
    if (!is.matrix(y) & !is(y, "dgCMatrix"))
        y <- .extract_counts(y, assay_sp, slot_sp)
    
    if (is.matrix(x)) {
        if (pnmf == "RcppML")
            # Make sure matrix is sparse
            x <- Matrix(x, sparse = TRUE)
    } else {
        if (pnmf == "NMF") 
            # Make sure matrix is dense
            x <- as.matrix(x)
    }
    
    # select genes in mgs or hvg
    if (!is.null(hvg)) {
        # Select union of genes between markers and HVG
        mod_genes <- union(unique(mgs[, gene_id]), hvg)
    } else {
        # Select genes from the marker genes only
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
    
    # Get seeding matrices
    if (verbose) message("Seeding initial matrices...")
    hw <- .init_nmf(x, groups, mgs, n_top, gene_id, group_id, weight_id)
    
    if (pnmf == "NMF") {
        .test_installed("NMF")
        # Seed NMF model
        seed <- NMF::nmfModel(W = hw$W, H = hw$H, model = paste0("NMF", model))
        # train NMF model
        if (verbose) message("Training NMF model...")
        mod <- NMF::nmf(x, rank, paste0(model, "NMF"), seed, ...)
    } else if (pnmf == "RcppML") {
        .test_installed("RcppML")
        if (verbose) message("Training NMF model...") 
        
        mod <- RcppML::nmf(
            data = x,
            k = rank,
            tol = tol,
            verbose = verbose,
            L1 = L1_nmf,
            L2 = L2_nmf,
            seed = hw$W)
        
        # Change nmfX to topic_X
        colnames(mod@w) <- paste0("topic_", seq_len(ncol(mod@w)))
        rownames(mod@h) <- paste0("topic_", seq_len(nrow(mod@h)))
    }
    # capture stop time
    t1 <- Sys.time()
    
    # print runtimes
    if (verbose) {
        dt <- round(difftime(t1, t0, units = "mins"), 2)
        message("Time for training: ", dt, "min")
    }
    
    # Extract NMFfit to list for consistency with RcppML
    mod <- .extract_nmf(mod, hw$W)
    
    # get topic profiles per cell type
    topic <- .topic_profiles(mod, groups)
    
    return(list("mod" = mod, "topic" = topic))
}
