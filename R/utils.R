#' @importFrom sparseMatrixStats rowSds
#' @importFrom Matrix t
.scale_uv <- function(x) {
    sds <- rowSds(x, na.rm = TRUE)
    # TODO find a more efficient way of scaling the matrix
    # t1 <- t(scale(t(x), center = FALSE, scale = sds))
    # Scale by gene (each row by its sd) for unit variance
    t1 <- x / sds
    t1
}

#' @importFrom Matrix Matrix
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
    # mgs <- lapply(ks, function(k) {
    #     g1 <- mgs[[k]][[gene_id]]
    #     g2 <- unlist(lapply(mgs[ks != k], `[[`, gene_id))
    #     mgs[[k]][!g1 %in% g2, , drop = FALSE]
    # })

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

#' .init_nmf <- function(x,
#'     groups,
#'     mgs,
#'     n_top = NULL,
#'     gene_id = "gene",
#'     group_id = "cluster",
#'     weight_id = "weight") {
#'     # check validity of input arguments
#'     if (is.null(n_top)) {
#'         n_top <- max(table(mgs[[group_id]]))
#'     }
#'     stopifnot(
#'         is.character(gene_id), length(gene_id) == 1,
#'         is.character(group_id), length(group_id) == 1,
#'         is.character(weight_id), length(weight_id) == 1,
#'         c(gene_id, group_id, weight_id) %in% names(mgs),
#'         is.numeric(n_top), length(n_top) == 1, round(n_top) == n_top)
#'     
#'     ng <- nrow(x)
#'     nc <- ncol(x)
#'     names(ks) <- ks <- unique(groups)
#'     
#'     # subset 'n_top' features
#'     mgs <- split(mgs, mgs[[group_id]])
#'     mgs <- lapply(mgs, function(df) {
#'         o <- order(df[[weight_id]], decreasing = TRUE)
#'         n <- ifelse(nrow(df) < n_top, nrow(df), n_top)
#'         df[o, ][seq_len(n), ]
#'     })
#'     
#'     # subset unique features
#'     mgs <- lapply(ks, function(k) {
#'         g1 <- mgs[[k]][[gene_id]]
#'         g2 <- unlist(lapply(mgs[ks != k], `[[`, gene_id))
#'         mgs[[k]][!g1 %in% g2, , drop = FALSE]
#'     })
#'     
#'     # W is of dimension (#groups)x(#features) with W(i,j)
#'     # equal to weight if j is marker for i, and ~0 otherwise
#'     W <- vapply(ks, function(k) {
#'         w <- numeric(ng) + 1e-12
#'         names(w) <- rownames(x)
#'         ws <- mgs[[k]][[weight_id]]
#'         w[mgs[[k]][[gene_id]]] <- ws
#'         return(w)
#'     }, numeric(ng))
#'     
#'     # there is no need to initialize H
#'     tp <- paste0("topic_", seq_len(length(ks)))
#'     dimnames(W) <- list(rownames(x), tp)
#'     return(W)
#' }

#' @importFrom Matrix Matrix rowSums
.filter <- function(x, y) {
    # remove undetected features
    .fil <- function(.) {
        i <- rowSums(.) > 0
        .[i, , drop = FALSE]
    }
    x <- .fil(x)
    
    # keep only shared features
    if (!is.null(y))
        x <- x[intersect(rownames(x), y), ]
    
    if (nrow(x) < 10) {
        stop(
            "Insufficient number of features shared",
            " between single-cell and mixture dataset.")
    }
    return(x)
}


#' @importFrom sparseMatrixStats colMedians
.topic_profiles <- function(mod, groups) {
    # Treat mod differently if it comes from NMF or RcppML
    df <- data.frame(t(mod$h))
    dfs <- split(df, groups)
    res <- vapply(
        dfs, function(df)
            colMedians(as.matrix(df)),
        numeric(ncol(df))
    )
    rownames(res) <- paste0("topic_", seq_len(nrow(res)))
    return(t(res))
}


#' @importFrom sparseMatrixStats rowSums2
.pred_hp <- function(
        x, mod, scale = TRUE, verbose = TRUE,
        L1_nnls = 0, L2_nnls = 0, threads = 0
    ) {
    W <- mod$w
    # remove all genes that are all 0s
    g0 <- rowSums2(x) > 0
    # Return a warning about genes being removed
    if (!all(g0) & verbose)
        message("Removing genes in mixture matrix that are all 0s")
    x <- x[g0, ]
    
    # Subset to shared genes between SP and SC
    if (verbose)
        message("Keep intersection of genes between W and mixture matrix")
    gi <- intersect(rownames(W), rownames(x))
    x <- x[gi, ]
    W <- W[gi, ]
    
    # Check there are enough shared features
    if (nrow(x) < 10) {
        stop(
            "Insufficient number of features, <10, shared",
            " between trained model and mixture dataset.")
    }
    if (scale) {
        x <- .scale_uv(x)
    }
    
    # TODO sometimes this can predict all to 0 if not scaled
    # If I do this we get the same since colSums(W) = 1 for all coummns
    # Use a very very mild regularization at this step
    # TODO revert back to native RCPP code works
    y <- predict_nmf(as(x, "dgCMatrix"), t(W), L1_nnls, L2_nnls, threads)
    # y <- RcppML::project(
    #   A = as(x, "dgCMatrix"),
    #   w = W,
    #   L1 = L1_nnls,
    #   nonneg = TRUE)
    
    # TODO set up a test to deal when a column in y is all 0s, meaning all the topics are 0 for that cell type
    
    # Assign names
    rownames(y) <- rownames(mod$h)
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

# Helper function to substitute the S4 method.
# This function takes in an object of class accepted in SPOTlight, it
# extracts the count/expression matrix specified and returns a matrix
.extract_counts <- function(x, assay, slot) {
    # Iterate over all the accepted classes and return expression matrix
    
    # Extract count matrix from object
    if (is(x, "Seurat")) {
        .test_installed(c("SeuratObject"))
        # Stop if there are no images or the name selected doesn't exist
        stopifnot(
            # Stop if there are no images
            !is.null(SeuratObject::Assays(x)),
            # Stop if the assay doesn't exist
            assay %in% SeuratObject::Assays(x)
        )
        
        # Extract Seurat coordinates
        x <- SeuratObject::GetAssayData(x, slot, assay)
    } else if (is(x, "SpatialExperiment") | is(x, "SingleCellExperiment")) {
        .test_installed(c("SummarizedExperiment"))
        
        # Stop if there are no images or the name selected doesn't exist
        stopifnot(
            # Stop if there are no images
            !is.null(SummarizedExperiment::assayNames(x)),
            # Stop if the image doesn't exist
            slot %in% SummarizedExperiment::assayNames(x),
            # Return error if there are no colnames in the object
            !is.null(colnames(x))
        )
        ## Extract SCE-SE coordinates
        x <- SummarizedExperiment::assay(x, slot)
    }
    
    # Process expression matrix
    if (is(x, "DelayedMatrix")) {
        # Convert to matrix
        rn <- rownames(x)
        cn <- colnames(x)
        x <- Matrix(x, sparse = TRUE, nrow = nrow(x), ncol = ncol(x))
        rownames(x) <- rn
        colnames(x) <- cn
    } else if (is(x, "dgCMatrix") | is.matrix(x)) {
        x
    } else {
        stop("Couldn't extract counts. Please check class(x) is a
        SingleCellExpriment, SpatialExperiment, Seurat, matrix, DelayedMatrix
        or dgCMatrix.")
    }
    return(x)
    
}

# Take an array representing an image and plot it with ggplot2
#' @import ggplot2
.plot_image <- function(x) {
    # Check necessary packages are installed and if not STOP
    .test_installed(c("grid", "ggplot2"))
    
    x <- grid::rasterGrob(x,
        interpolate = FALSE,
        width = grid::unit(1, "npc"),
        height = grid::unit(1, "npc"))
    
    ggplot() +
        annotation_custom(
            grob = x,
            xmin = 0,
            xmax = ncol(x$raster),
            ymin = 0,
            ymax = nrow(x$raster)) + 
        coord_fixed(
            xlim = c(0, ncol(x$raster)),
            ylim = c(0, nrow(x$raster))) + 
        theme_void()
        # theme_classic()
}

# Extract image and convert it to array from allowed classes
.extract_image <- function(x, slice = NULL) {
    # Iterate over all the accepted classes and convert the image to array
    if (is.character(x)) {
        .test_installed(c("jpeg", "png"))
        
        # Check if the file exists
        stopifnot(file.exists(x))
        
        # Check the file is in the right format
        typ <- c("jpg", "jpeg", "png")
        pat <- paste0(".", typ, "$")
        idx <- vapply(pat, grepl, x = x, logical(1), ignore.case = TRUE)
        if (!any(idx)) {
            stop("'x' should be of file type JPG, JPEG or PNG")
        }
        
        # Read file
        x <- switch(typ[idx],
            png = png::readPNG(x),
            jpeg::readJPEG(x))
        
    } else if (is(x, "Seurat")) {
        .test_installed(c("SeuratObject"))
        # Stop if there are no images or the name selected doesn't exist
        stopifnot(
            !is.null(SeuratObject::Images(x)),
            slice %in% SeuratObject::Images(x))
        
        # If image is null use the first slice
        if (is.null(slice)) 
            slice <- SeuratObject::Images(x)[1]
        
        # Extract Image in raster format
        x <- SeuratObject::GetImage(x, image = slice, mode = "raster")
        # Conver to matrix
        x <- as.matrix(x)
        
    } else if (is(x, "SpatialExperiment")) {
        
        .test_installed(c("SpatialExperiment"))
        
        # Stop if there are no images or the name selected doesn't exist
        stopifnot(
            !is.null(SpatialExperiment::getImg(x)),
            slice %in% SpatialExperiment::imgData(x)[1, "sample_id"]
        )
        
        # If image is null use the first slice
        if (is.null(slice)) 
            slice <- SpatialExperiment::imgData(x)[1, "sample_id"]
        
        # Convert to raster
        x <- SpatialExperiment::imgRaster(x, sample_id = slice)
        x <- as.matrix(x)
    } else {
        stop("Couldn't extract image, See ?plotImage for valid image inputs.")
    }
    return(x)
}

# When assigning cells to groups in trainNMF and SPOTlight if groups is set to
# NULL use the cell identities/labels. If it is not a Seurat or SCE return error
.set_groups_if_null <- function(x) {
    ## Seurat ##
    if (is(x, "Seurat")) {
        # Extract idents
        idents <- SeuratObject::Idents(x)
        if (is.null(idents)) {
            stop("SeuratObject::Idents(x) is NULL")
        } else {
            warning("Grouping cells into celltypes by Idents(x)")
            groups <- as.character(idents)
        }
        
    ## SCE ##
    } else if (is(x, "SingleCellExperiment")) {
        # Extract idents
        idents <- SingleCellExperiment::colLabels(x)
        if (is.null(idents)) {
            stop("SingleCellExperiment::colLabels(x) is NULL")
        } else {
            warning("Grouping cells into celltypes
                    by SingleCellExperiment::colLabels(x)")
            groups <- as.character(idents)
        }
    ## other ##
    } else {
        stop("Parameter groups needs to be defined.")
    }
    groups
}

# Helper function to extract elements of interest from objects NMFfit
# (NMF package) and nmf (RcppML) and returns a list with relevant information
# consistent between both of them
# .extract_nmf <- function(mod, smtx) {
#     if (is(mod, "NMFfit")) {
#         mod <- list(
#             "w" = NMF::basis(mod),
#             "d" = NULL,
#             "h" = NMF::coef(mod),
#             "misc" = list(
#                 "tol" = NULL,
#                 "iter" = mod@extra$iteration,
#                 "runtime" = mod@runtime,
#                 "mse" = NULL,
#                 "w_init" = smtx)
#         )
#     } else if (is.list(mod)) {
#         mod <- list(
#             "w" = mod$w,
#             "d" = mod$d,
#             "h" = mod$h,
#             "misc" = list(
#                 "tol" = NULL,
#                 "iter" = NULL,
#                 "runtime" = NULL,
#                 "mse" = NULL,
#                 "w_init" = NULL)
#         )
#     } else {
#         stop("mod is neither an 'NMFfit' or 'nmf' object ")
#     }
#     
#     return(mod)
# }
