#' @rdname data
#' @name data
#' @aliases mockSC mockSP getMGS
#' @title Synthetic single-cell, mixture and marker data
#'
#' @description
#' \code{mockSC/mockSP()} are designed to generate synthetic single-cell and
#' spatial mixture data. These data are not meant to represent biologically
#' meaningful use-cases, but are solely intended for use in examples, for
#' unit-testing, and to demonstrate \code{SPOTlight}'s general functionality.
#' Finally, \code{.get_mgs()} implements a statistically naive way to select
#' markers from single-cell data; again, please don't use it in real life.
#'
#' @param ng,nc,nt,ns integer scalar specifying the number
#'   of genes, cells, types (groups) and spots to simulate.
#' @param n_top integer scalar specifying the number
#'   of markers to select per group.
#'
#' @return
#' \itemize{
#' \item{\code{mockSC} returns a \code{SingleCellExperiment}
#'   with rows = genes, columns = single cells, and cell metadata
#'   (\code{colData}) column \code{type} containing group identifiers.}
#' \item{\code{mockSP} returns a \code{SingleCellExperiment}
#'   with rows = genes, columns = single cells, and cell metadata
#'   (\code{colData}) column \code{type} containing group identifiers.}
#' \item{\code{getMGS} returns a \code{data.frame} with \code{nt*n_top}
#'   rows and 3 columns: gene and type (group) identifier, as well as the
#'   gene's weight = the proportion of counts accounted for by that type.}
#' }
#'
#' @examples
#' sce <- mockSC()
#' spe <- mockSP(sce)
#' mgs <- getMGS(sce)
NULL

#' @rdname data
#' @importFrom SingleCellExperiment cbind SingleCellExperiment
#' @importFrom stats rnbinom runif
#' @export
mockSC <- function(ng = 200, nc = 50, nt = 3) {
    z <- lapply(seq_len(nt), \(t) {
        ms <- 2^runif(ng, 2, 10)
        ds <- 0.5 + 100 / ms
        y <- rnbinom(ng * nc, mu = ms, size = 1 / ds)
        y <- matrix(y, nrow = ng, ncol = nc)
        dimnames(y) <- list(
            paste0("gene", seq_len(ng)),
            paste0("cell", seq_len(nc))
        )
        x <- SingleCellExperiment(list(counts = y))
        x$type <- factor(
            paste0("type", t),
            paste0("type", seq_len(nt))
        )
        return(x)
    })
    do.call(cbind, z)
}

#' @rdname data
#' @param x Single cell experiment object ç
#' @importFrom Matrix rowSums
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
mockSP <- function(x, ns = 100) {
    z <- replicate(ns, {
        # sample number of cells
        nc <- sample(5, 1)
        # sample reference cells
        cs <- sample(ncol(x), nc)
        # sum up counts & rescale
        y <- counts(x[, cs])
        y <- rowSums(y)
        # compute composition
        n <- table(x$type[cs]) / nc
        n <- c(unclass(n))
        list(y, n)
    })
    # get counts
    y <- t(do.call(rbind, z[1, ]))
    dimnames(y) <- list(
        rownames(x),
        paste0("spot", seq_len(ns))
    )
    # get compositions
    fq <- do.call(rbind, z[2, ])
    rownames(fq) <- colnames(y)
    # sample coordinates
    xy <- matrix(runif(2 * ns), ncol = 2)
    dimnames(xy) <- list(colnames(y), c("x", "y"))
    SingleCellExperiment(
        list(counts = y),
        colData = data.frame(xy),
        metadata = list(props = fq)
    )
}

#' @rdname data
#' @param n_top integer specifying the number of  
#'   marker genes to extract for each cluster.
#' @importFrom Matrix colSums rowSums
#' @importFrom SingleCellExperiment counts
#' @importFrom stats aggregate
#' @export
getMGS <- function(x, n_top = 10) {
    # compute sum of counts by group
    y <- aggregate(t(counts(x)), list(x$type), sum)
    rownames(y) <- y[, 1]
    # Remove group column
    y <- t(y[, -1])
    # get proportion of counts by group
    z <- lapply(rownames(y), \(gene) {
        p <- prop.table(y[gene, ])
        i <- which.max(p)
        type <- names(i)
        weight <- p[i]
        data.frame(gene, type, weight)
    })
    z <- do.call(rbind, z)
    rownames(z) <- NULL
    # select 'top_n' in each group
    z <- split(z, z$type)
    # Iterate over groups and sort within them
    z <- lapply(z, \(.) {
        # Get indexes of the positions in the sorted order
        o <- order(.$weight, decreasing = TRUE)
        # order the markers
        n <- nrow(.)
        if (n < n_top)
            n_top <- n
        .[o, ][seq_len(n_top), ]
    })
    z <- do.call(rbind, z)
    rownames(z) <- NULL
    return(z)
}
