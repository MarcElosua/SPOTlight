#' @importFrom matrixStats rowSds
.scale_uv <- function(x)
{
    sds <- rowSds(x, na.rm = TRUE)
    t(scale(t(x), center = FALSE, scale = sds))
}

.init_nmf <- function(x, 
    groups,
    mgs, 
    n_top = NULL, 
    gene_id = "gene", 
    group_id = "cluster", 
    weight_id = "weight")
{
    # check validity of input arguments
    if (is.null(n_top)) 
        n_top <- max(table(mgs[[group_id]]))
    stopifnot(
        is.character(gene_id), length(gene_id) == 1,
        is.character(group_id), length(group_id) == 1,
        is.character(weight_id), length(weight_id) == 1,
        c(gene_id, group_id, weight_id) %in% names(mgs),
        is.numeric(n_top), length(n_top) == 1, round(n_top) == n_top)

    ng <- nrow(x)
    nc <- ncol(x)
    names(ks) <- ks <- unique(groups)

    # subset 'n_top' features
    mgs <- split(mgs, mgs[[group_id]])
    mgs <- lapply(mgs, \(df) {
        o <- order(df[[weight_id]], decreasing = TRUE)
        n <- ifelse(nrow(df) < n_top, nrow(df), n_top)
        df[o, ][seq_len(n), ]
    })
    
    # subset unique features
    mgs <- lapply(ks, \(k) {
        g1 <- mgs[[k]][[gene_id]]
        g2 <- unlist(lapply(mgs[ks != k], `[[`, gene_id))
        mgs[[k]][!g1 %in% g2, , drop = FALSE]
    })
    
    # W is of dimension (#groups)x(#features) with W(i,j)
    # equal to weight if j is marker for i, and ~0 otherwise
    W <- vapply(ks, \(k) {
        w <- numeric(ng) + 1e-12
        names(w) <- rownames(x)
        ws <- mgs[[k]][[weight_id]]
        w[mgs[[k]][[gene_id]]] <- ws
        return(w)
    }, numeric(ng))
    
    # H is of dimension (#groups)x(#samples) with H(i,j)
    # equal to 1 if j is in i, and ~0 otherwise
    cs <- split(seq_len(nc), groups)
    H <- t(vapply(ks, \(k) {
        h <- numeric(nc) + 1e-12
        h[cs[[k]]] <- 1
        return(h)
    }, numeric(nc)))
    
    dimnames(W) <- list(rownames(x), ks)
    dimnames(H) <- list(ks, colnames(x))
    return(list("W" = W, "H" = H))
}

.filter <- function(x, y) 
{
    # remove undetected features
    .fil <- \(.) {
        i <- rowSums(.) > 0
        .[i, , drop = FALSE]
    }
    x <- .fil(x)
    y <- .fil(y)
    
    # keep only shared features
    i <- intersect(
        rownames(x), 
        rownames(y))
    if (length(i) < 10)
        stop("Insufficient number of features shared",
            " between single-cell and mixture dataset.")
    return(x[i, ])
}

#' @importFrom Matrix rowSums
#' @importFrom NMF nmf nmfModel
.train_nmf <- function(x, y,
    groups,
    mgs,
    n_top = NULL,
    gene_id = "gene",
    group_id = "cluster",
    weight_id = "weight",
    model = c("ns", "std"),
    scale = TRUE,
    verbose = TRUE)
{
    # check validity of input arguments
    model <- match.arg(model)
    
    # drop features that are undetected 
    # in single-cell and/or mixture data
    x <- .filter(x, y)
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
    # TODO: figure out if these models can be
    # different / if it matters if they are
    
    # train NMF model
    if (verbose) message("Training NMF model")
    mod <- nmf(x, rank, paste0(model, "NMF"), seed) 

    # capture top time
    t1 <- Sys.time()
    
    # print runtimes
    if (verbose) {
        dt <- round(difftime(t1, t0, units = "mins"), 2)
        message("Time for training: ", dt, "min")
    }
    return(mod)
}

#' @importFrom matrixStats colMedians
#' @importFrom NMF coef
.topic_profiles <- function(mod, groups)
{
    df <- data.frame(t(coef(mod)))
    dfs <- split(df, groups)
    res <- vapply(dfs, \(df) 
        colMedians(as.matrix(df)),
        numeric(ncol(df)))
    rownames(res) <- names(dfs)
    return(t(res))
}

#' @importFrom NMF basis
#' @importFrom nnls nnls
.pred_prop <- function(x, mod, scale = TRUE, verbose = TRUE)
{
    # TODO: 
    # if 'scale = TRUE' in 'SPOTlight()', this is already 
    # done by '.train_nmf()'. could be removed here?
    W <- basis(mod)
    x <- x[rownames(W), ]
    if (scale) 
        x <- .scale_uv(x)
    y <- vapply(
        seq_len(ncol(x)), 
        \(i) nnls(W, x[, i])$x,
        numeric(ncol(W)))
    rownames(y) <- dimnames(mod)[[3]]
    colnames(y) <- colnames(x)
    return(y)
}

#' @importFrom nnls nnls
.deconvolute <- function(x, mod, ref, scale = TRUE, min_prop = 0.01, verbose = TRUE) 
{
    mat <- .pred_prop(x, mod, scale)
    if (verbose) message("Deconvoluting mixture data")
    res <- vapply(seq_len(ncol(mat)), \(i) {
        pred <- nnls(ref, mat[, i])
        prop <- prop.table(pred$x)
        # drop groups that fall below 'min_prop' & update
        prop[prop < min_prop] <- 0
        prop <- prop.table(prop)
        # compute residual sum of squares
        ss <- sum(mat[, i]^2)
        # compute percentage of unexplained residuals
        err <- pred$deviance/ss
        c(prop, err)
    }, numeric(ncol(ref)+1))
    # set dimension names
    rownames(res) <- c(dimnames(mod)[[3]], "res_ss")
    colnames(res) <- colnames(mat)
    return(t(res))
}
