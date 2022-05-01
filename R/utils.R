#' @importFrom sparseMatrixStats rowSds
#' @importFrom Matrix t
.scale_uv <- function(x) {
    sds <- rowSds(x, na.rm = TRUE)
    # TODO find a more efficient way of scaling the matrix
    # t1 <- t(scale(t(x), center = FALSE, scale = sds))
    t1 <- t(t(x) / sds)
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

# Helper function to substitute the S4 method.
# This function takes in an object of class accepted in SPOTlight, it
# extracts the count/expression matrix specified and returns a matrix
.extract_counts <- function(x, assay, slot) {
    # Iterate over all the accepted classes and return expression matrix
    if (is(x, "DelayedMatrix")) {
        # Convert to matrix
        rn <- rownames(x)
        cn <- colnames(x)
        x <- Matrix(x, sparse = TRUE, nrow = nrow(x), ncol = ncol(x))
        rownames(x) <- rn
        colnames(x) <- cn
    } else if (is(x, "Seurat")) {
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
