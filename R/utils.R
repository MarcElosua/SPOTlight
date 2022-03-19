#' @importFrom sparseMatrixStats rowSds
#' @importFrom Matrix t
.scale_uv <- function(x) {
    sds <- rowSds(x, na.rm = TRUE)
    # TODO find a more efficient way of scaling the matrix
    # t1 <- t(scale(t(x), center = FALSE, scale = sds))
    t1 <- t(t(x) / sds)
    print(t1[1:5, 1:5])
    print(is(t1))
    t1
}

.init_nmf <- function(x,
    groups,
    mgs,
    n_top = NULL,
    gene_id = "gene",
    group_id = "cluster",
    weight_id = "weight") {
    # check validity of input arguments
    if (is.null(n_top)) {
        n_top <- max(table(mgs[[group_id]]))
    }
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
    mgs <- lapply(mgs, function(df) {
        o <- order(df[[weight_id]], decreasing = TRUE)
        n <- ifelse(nrow(df) < n_top, nrow(df), n_top)
        df[o, ][seq_len(n), ]
    })

    # subset unique features
    mgs <- lapply(ks, function(k) {
        g1 <- mgs[[k]][[gene_id]]
        g2 <- unlist(lapply(mgs[ks != k], `[[`, gene_id))
        mgs[[k]][!g1 %in% g2, , drop = FALSE]
    })

    # W is of dimension (#groups)x(#features) with W(i,j)
    # equal to weight if j is marker for i, and ~0 otherwise
    W <- vapply(ks, function(k) {
        w <- numeric(ng) + 1e-12
        names(w) <- rownames(x)
        ws <- mgs[[k]][[weight_id]]
        w[mgs[[k]][[gene_id]]] <- ws
        return(w)
    }, numeric(ng))

    # H is of dimension (#groups)x(#samples) with H(i,j)
    # equal to 1 if j is in i, and ~0 otherwise
    cs <- split(seq_len(nc), groups)
    H <- t(vapply(ks, function(k) {
        h <- numeric(nc) + 1e-12
        h[cs[[k]]] <- 1
        return(h)
    }, numeric(nc)))
    
    tp <- paste0("topic_", seq_len(length(ks)))
    dimnames(W) <- list(rownames(x), tp)
    dimnames(H) <- list(tp, colnames(x))
    return(list("W" = W, "H" = H))
}

.filter <- function(x, y) {
    # remove undetected features
    .fil <- function(.) {
        i <- rowSums(.) > 0
        .[i, , drop = FALSE]
    }
    x <- .fil(x)
    y <- .fil(y)

    # keep only shared features
    i <- intersect(
        rownames(x),
        rownames(y))
    
    if (length(i) < 10) {
        stop(
            "Insufficient number of features shared",
            " between single-cell and mixture dataset.")
    }
    return(x[i, ])
}

#' @importFrom sparseMatrixStats colMedians
#' @importFrom NMF coef
.topic_profiles <- function(mod, groups) {
    # Treat mod differently if it comes from NMF or RcppML
    if (is(mod, "NMFfit")) {
        df <- data.frame(t(coef(mod)))
    } else if (is.list(mod)) {
        df <- data.frame(t(mod$h))
    }
    
    dfs <- split(df, groups)
    res <- vapply(
        dfs, function(df)
            colMedians(as.matrix(df)),
        numeric(ncol(df))
    )
    rownames(res) <- paste0("topic_", seq_len(nrow(res)))
    return(t(res))
}

#' @importFrom NMF basis
#' @importFrom nnls nnls
.pred_prop <- function(x, mod, scale = TRUE, verbose = TRUE) {
    # Keep basis sparse
    if (class(mod, "NMFfit")) {
        W <- data.frame(basis(mod))
    } else if (is.list(mod)) {
        W <- data.frame(mod$w)
    }
    
    x <- x[rownames(W), ]
    if (scale) {
        x <- .scale_uv(x)
    }
    # TODO go from here
    # Error in nnls(W, x[, i]) : 
    # 'list' object cannot be coerced to type 'double'
    y <- vapply(
        seq_len(ncol(x)), 
        function(i) nnls(W, x[, i])$x,
        numeric(ncol(W)))
    
    rownames(y) <- dimnames(mod)[[3]]
    colnames(y) <- colnames(x)
    return(y)
}

# Test if a package is installed
# x is a stringr or vector of strings of packages names
# to test if they are installed
.test_installed <- function(x) {
    # Check which packages aren't installed
    t <- vapply(x, function(i)
        isFALSE(requireNamespace(i, quietly = TRUE)), numeric(1))
    x <- x[t == 1]

    if (length(x) > 0) {
        x <- paste(x, collapse = ", ")
        stop("Please install package/s: ", x)
    }
}
