#' @name runDeconvolution
#' @rdname runDeconvolution
#' @title Run Deconvolution using NNLS model
#'
#' @aliases runDeconvolution
#'
#' @description This function takes in the mixture data, the trained model & the
#'   topic profiles and returns the proportion of each cell type within each
#'    mixture
#'
#' @param x mixture dataset. Can be a numeric matrix,
#'   \code{SingleCellExperiment}, \code{SpatialExperiment} or
#'   \code{SeuratObjecy}.

#' @param mod object as obtained from trainNMF.
#' @param ref object of class matrix containing the topic profiles for each cell
#'  type as obtained from trainNMF.
#' @param assay if the object is of Class \code{Seurat}, character string
#'   specifying the assay from which to extract the expression matrix.
#'   By default "Spatial".
#' @param slot if the object is of Class \code{Seurat}, character string
#'   specifying the slot from which to extract the expression matrix. If the
#'   object is of class \code{SpatialExperiment} indicates matrix to use.
#'   By default "counts".
#' @inheritParams SPOTlight
#'
#' @return base a list where the first element is a list giving the NMF model and
#'   the second is a matrix containing the topic profiles learnt.
#'
#' @author Marc Elosua Bayes, Zach DeBruine, and Helena L Crowell
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
#' # Run deconvolution
#' decon <- runDeconvolution(
#'     x = spe,
#'     mod = res[["mod"]],
#'     ref = res[["topic"]])
NULL

#' @rdname runDeconvolution
#' @export
runDeconvolution <- function(
    x,
    mod,
    ref,
    scale = TRUE,
    min_prop = 0.01,
    verbose = TRUE,
    assay = "Spatial",
    slot = "counts",
    L1_nnls_topics = 0,
    L2_nnls_topics = 0,
    L1_nnls_prop = 0,
    L2_nnls_prop = 0,
    threads = 0,
    ...) {

    # Class checks
    stopifnot(
        # Check x inputs
        is.matrix(x) | is(x, "DelayedMatrix") | is(x, "dgCMatrix") |
            is(x, "Seurat") | is(x, "SingleCellExperiment") |
            is(x, "SpatialExperiment"),
        # Check mod inputs
        is.list(mod),
        # check ref
        is.matrix(ref),
        # Check assay name
        is.character(assay), length(assay) == 1,
        # Check slot name
        is.character(slot), length(slot) == 1,
        # Check scale and verbose
        is.logical(scale), length(scale) == 1,
        is.logical(verbose), length(verbose) == 1,
        # Check min_prop numeric
        is.numeric(min_prop), length(min_prop) == 1,
        min_prop >= 0, min_prop <= 1
    )

    # Extract expression matrix
    if (!is.matrix(x))
        x <- .extract_counts(x, assay, slot)

    # Get topic profiles for mixtures
    mat <- .pred_hp(
        x = x, mod = mod, scale = scale, verbose = verbose,
        L1_nnls = L1_nnls_topics, L2_nnls = L2_nnls_topics, threads = threads)
    
    if (verbose) message("Deconvoluting mixture data...")
    ref_scale <- t(t(ref) / colSums(ref))
    # ref_scale <- t(ref) / colSums(ref)
    
    # TODO I want the below line to do but it doesn't!
    # pred <- t(mat) %*% t(ref_scale)
    # res <- prop.table(pred, 1)
    # TODO come back to change this with the native RCPP code
    # pred <- predict_nmf(
    #     A_ = as(mat, "dgCMatrix"),
    #     w = ref_scale,
    #     L1 = L1_nnls,
    #     L2 = L2_nnls,
    #     threads = threads)
    pred <- RcppML::project(
      A = as(mat, "dgCMatrix"),
      w = t(ref_scale),
      L1 = 0,
      nonneg = TRUE)
    # rownames(pred) <- rownames(ref_scale)
    # colnames(pred) <- colnames(mat)
    rownames(pred) <- rownames(ref_scale)
    colnames(pred) <- colnames(mat)

    # Proportions within each spot
    res <- prop.table(pred, 2)

    # TODO Check computation is correct for residuals
    # 1- t(ref_scale) %*% pred map pred to mat using ref_scale
    # 2- Check the differences between the original and re-mapped matrix
    # 3- sum the errors for each spot (column)
    err_mat <- (mat - t(ref_scale) %*% pred)^2
    err <- colSums(err_mat)
    names(err) <- colnames(res)

    return(list("mat" = t(res), "res_ss" = err))
}

