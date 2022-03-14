#' @importFrom nnls nnls
runDeconvolution <- function(x, mod, ref, scale = TRUE,
    min_prop = 0.01, verbose = TRUE) {
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
}
