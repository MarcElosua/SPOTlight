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
#' @param x single-cell dataset. Can be a numeric matrix, Can be a
#'   numeric matrix, \code{SingleCellExperiment} or \code{SeuratObjecy}.
#' @param y Null if you want to train the model with all the genes in the SC
#'    data or a character vector with the rownames of the mixture dataset to 
#'    subset the gene set used to the intersection between them.
#' @param assay_sc if the object is of Class \code{Seurat}, character string
#'   specifying the assay from which to extract the expression matrix.
#'   By default "RNA".
#' @param slot_sc if the object is of Class \code{Seurat}, character string
#'   specifying the slot from which to extract the expression matrix. If the
#'   object is of class \code{SingleCellExperiment} indicates matrix to use.
#'   By default "counts".

#' @inheritParams SPOTlight
#'
#' @return a list where the first element is a list with the NMF model
#'   information and the second is a matrix containing the topic profiles
#'   learnt per cell type.
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
#'     y = rownames(spe),
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

# Key here to load t & Matrix so sparse matrices can be transposed
#' @importFrom Matrix Matrix t
#' @export
trainNMF <- function(
    x,
    y = NULL,
    groups = NULL,
    mgs,
    n_top = NULL,
    gene_id = "gene",
    group_id = "cluster",
    weight_id = "weight",
    hvg = NULL,
    scale = TRUE,
    verbose = TRUE,
    L1_nmf = 0,
    L2_nmf = 0,
    tol = 1e-05,
    maxit = 100,
    threads = 0,
    assay_sc = "RNA",
    slot_sc = "counts",
    ...) {
    
    if (is.null(n_top))
        n_top <- max(table(mgs[[group_id]]))
    ids <- c(gene_id, group_id, weight_id)
    
    # convert mgs to dataframe if it is not already
    if (!is.data.frame(mgs)) {
      # check.names=FALSE to ensure the ids specified by the user are unchanged
      mgs <- data.frame(mgs, check.names = FALSE)
    }
    
    stopifnot(
        is.numeric(x) | is(x, "dgCMatrix") |
            is(x, "Seurat") | is(x, "SingleCellExperiment") |
            is(x, "DelayedMatrix"), 
        (is.vector(y) & is.character(y)) | is.null(y),
        is.character(ids), length(ids) == 3, ids %in% names(mgs),
        is.null(groups) | length(groups) == ncol(x),
        is.logical(scale), length(scale) == 1,
        is.logical(verbose), length(verbose) == 1,
        is.numeric(L1_nmf), length(L1_nmf) == 1,
        is.numeric(L2_nmf), length(L2_nmf) == 1,
        is.numeric(tol), length(tol) == 1)
    
    # Set groups if x is SCE or SE and groups is NULL 
    if (is.null(groups))
        groups <- .set_groups_if_null(x)
    
    groups <- as.character(groups)

    # Check mgs is a dataframe or conver it to a df
    if (!is.data.frame(mgs)) {
        if (is(mgs, "tibble") || is(mgs, "list")) {
            mgs <- as.data.frame(mgs)
        } else stop("'mgs' should be a 'data.frame'")
    }
    
    # Stop if at least one of the groups doesn't have marker genes
    stopifnot(groups %in% mgs[[group_id]])
    
    # Extract expression matrices for x and y
    if (!is.matrix(x) & !is(x, "dgCMatrix"))
        x <- .extract_counts(x, assay_sc, slot_sc)
    
    # Make sure matrix is sparse
    # convert matrix to dgCMatrix, 
    # if it is already then nothing is done
    x <- as(x, "dgCMatrix")
    
    # Set y no rownames X if NULL
    if (is.null(y))
        y <- rownames(x)
    
    # select genes in mgs or hvg
    if (!is.null(hvg)) {
        # Select union of genes between markers and HVG
        mod_genes <- union(unique(mgs[[gene_id]]), hvg)
    } else {
        # Select genes from the marker genes only
        mod_genes <- unique(mgs[[gene_id]])
    }
    
    # Select intersection between interest and present in x (sce) & y (spe)
    mod_genes <- intersect(mod_genes, intersect(rownames(x), y))
    
    # drop features that are undetected in single-cell and/or mixture data
    x <- .filter(x[mod_genes, ], y)
    
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
    if (verbose) message("Seeding NMF model...")
    hw <- .init_nmf(x, groups, mgs, n_top, gene_id, group_id, weight_id)
    # w_init <- .init_nmf(x, groups, mgs, n_top, gene_id, group_id, weight_id)
    
    if (verbose) message("Training NMF model...") 
    
    # call to C++ routine
    mod <- run_nmf(x, t(x), tol, maxit, verbose, L1_nmf, L2_nmf, threads, t(hw$W))
    
    # Change nmfX to topic_X
    colnames(mod$w) <- paste0("topic_", seq_len(ncol(mod$w)))
    rownames(mod$h) <- paste0("topic_", seq_len(nrow(mod$h)))
    rownames(mod$w) <- rownames(x)
    colnames(mod$h) <- colnames(x)

    
    # capture stop time
    t1 <- Sys.time()
    
    # print runtimes
    if (verbose) {
        dt <- round(difftime(t1, t0, units = "mins"), 2)
        message("Time for training: ", dt, "min")
    }
    
    # Extract NMFfit to list for consistency with RcppML
    # mod <- .extract_nmf(mod, hw$W)
    
    # get topic profiles per cell type
    topic <- .topic_profiles(mod, groups)
    
    return(list("mod" = mod, "topic" = topic))
}
