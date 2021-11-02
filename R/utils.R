#' @importFrom methods is
#' @importFrom scuttle logNormCounts
#' @importFrom scran getTopHVGs modelGeneVar
#' @importFrom SummarizedExperiment rowData
.downsample_sce <- function(sce, n_cells = 100, n_genes = 3e3) {
    # check validity of input arguments
    stopifnot(
        is(sce, "SingleCellExperiment"),
        is.numeric(n_cells), 
        length(n_cells) == 1, 
        round(n_cells) == n_cells,
        is.numeric(n_genes), 
        length(n_genes) == 1, 
        round(n_genes) == n_genes)
    
    if (is.numeric(n_genes)) {
        stopifnot(
            n_genes > 0, round(n_genes) == n_genes, 
            length(n_genes) %in% c(1, ncol(sce)))
        if (length(n_genes) == 1) {
            # single numeric
            sce <- logNormCounts(sce)
            mat <- modelGeneVar(sce)
            gs <- getTopHVGs(mat, n = n_genes)
        } else {
            # numeric vector
            stopifnot(max(n_genes) <= nrow(sce))
        }
    } else if (is.logical(n_genes)) {
        # logical vector
        stopifnot(
            !isFALSE(n_genes), 
            length(n_genes) %in% c(1, nrow(sce)))
    } else if (is.character(n_genes)) {
        # metadata column
        gs <- rowData(sce)$n_genes
        stopifnot(!is.null(gs))
    }
    
    cs <- seq_len(ncol(sce))
    cs <- split(cs, colLabels(sce))
    cs <- lapply(cs, \(.) {
        n <- length(.)
        if (n_cells < n)
            n <- n_cells
        sample(., n)
    })
    cs <- unlist(cs)
    
    return(sce[gs, cs])
}
    
#' @importFrom Matrix rowSums
#' @importFrom matrixStats rowSds
#' @importFrom NMF nmf nmfModel
#' @importFrom SingleCellExperiment colLabels
#' @importFrom SummarizedExperiment assay
.train_nmf <- function(
    sce,
    spe,
    mgs,
    n_top = NULL,
    gene_id = "gene",
    group_id = "cluster",
    weight_id = "mean.AUC",
    model = c("ns", "std"),
    scale = TRUE,
    verbose = TRUE)
{
    # check validity of input arguments
    model <- match.arg(model)
    
    # remove undetected features
    .fil <- \(x) {
        y <- assay(x)
        i <- rowSums(y) > 0
        x[i, , drop = FALSE]
    }
    sce <- .fil(sce)
    spe <- .fil(spe)
    
    # keep only shared features
    gs <- intersect(
        rownames(sce), 
        rownames(spe))
    if (length(gs) < 10)
        stop("Insufficient number of features shared",
            " between single-cell and mixture dataset.")
    
    # subset single-cell data & markers
    sce <- sce[gs, ]
    mgs <- mgs[mgs[[gene_id]] %in% gs, ]
    
    # get counts & un-sparse them
    y <- assay(sce)
    if (!is.matrix(y))
        y <- as.matrix(y)
    
    # scale to unit variance (optional)
    if (scale) {
        if (verbose) 
            message("Scaling count matrix")
        sds <- rowSds(y, na.rm = TRUE)
        assay(sce) <- t(scale(t(y), center = FALSE, scale = sds))
    }
    
    # capture start time
    t0 <- Sys.time()
    
    # set model rank to number of groups
    rank <- length(unique(colLabels(sce)))
    
    # get seeding matrices (optional)
    seed <- if (TRUE) {
        if (verbose) 
            message("Seeding initial matrices")
        hw <- .init_nmf(sce, mgs, n_top, gene_id, group_id, weight_id)
        nmfModel(W = hw$W, H = hw$H, model = paste0("NMF", model))
    }
    # TODO: figure out if these models can be
    # different / if it matters if they are
    
    # train NMF model
    if (verbose)
        message("Training NMF model")
    mod <- nmf(y, rank, paste0(model, "NMF"), seed) 

    # capture top time
    t1 <- Sys.time()
    
    # print runtimes
    if (verbose) {
        dt <- round(difftime(t1, t0, units = "mins"), 2)
        message("Time for training: ", dt, "min")
    }

    return(mod)
}

.init_nmf <- function(sce, 
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
        is(sce, "SingleCellExperiment"),
        is.character(gene_id), length(gene_id) == 1,
        is.character(group_id), length(group_id) == 1,
        is.character(weight_id), length(weight_id) == 1,
        c(gene_id, group_id, weight_id) %in% names(mgs),
        is.numeric(n_top), length(n_top) == 1, round(n_top) == n_top)
    
    #y <- t(assay(sce))
    ng <- nrow(sce)
    nc <- ncol(sce)
    names(ks) <- ks <- unique(colLabels(sce))

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
        names(w) <- rownames(sce)
        ws <- mgs[[k]][[weight_id]]
        w[mgs[[k]][[gene_id]]] <- ws
        return(w)
    }, numeric(ng))
    
    # H is of dimension (#groups)x(#samples) with H(i,j)
    # equal to 1 if j is in i, and ~0 otherwise
    cs <- split(seq_len(nc), colLabels(sce))
    H <- t(vapply(ks, \(k) {
        h <- numeric(nc) + 1e-12
        h[cs[[k]]] <- 1
        return(h)
    }, numeric(nc)))
    
    dimnames(W) <- list(rownames(sce), ks)
    dimnames(H) <- list(ks, colnames(sce))
    return(list("W" = W, "H" = H))
}

#' @importFrom matrixStats colMedians
.topic_profiles <- function(sce, mod)
{
    df <- data.frame(t(coef(mod)))
    dfs <- split(df, colLabels(sce))
    res <- vapply(dfs, \(df) 
        colMedians(as.matrix(df)),
        numeric(ncol(df)))
    rownames(res) <- names(dfs)
    return(t(res))
}

.pred_prop <- function(mod, spe, scale)
{
    
}

.deconvolute <- function() 
{
    
}