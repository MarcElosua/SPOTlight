#' @rdname data
#' @name data
#' @aliases .mock_sc .mock_sp .get_mgs
#' @title Synthetic single-cell, mixture and marker data
#'
#' @description
#' \code{.mock_sc/sp()} are designed to generate synthetic single-cell and
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
#' \item{\code{.mock_sc} returns a \code{SingleCellExperiment}
#'   with rows = genes, columns = single cells, and cell metadata
#'   (\code{colData}) column \code{type} containing group identifiers.}
#' \item{\code{.mock_sp} returns a \code{SingleCellExperiment}
#'   with rows = genes, columns = single cells, and cell metadata
#'   (\code{colData}) column \code{type} containing group identifiers.}
#' \item{\code{.get_mgs} returns a \code{data.frame} with \code{nt*n_top}
#'   rows and 3 columns: gene and type (group) identifier, as well as the
#'   gene's weight = the proportion of counts accounted for by that type.}
#' }
#'
#' @examples
#' sce <- .mock_sc()
#' spe <- .mock_sp(sce)
#' mgs <- .get_mgs(sce)
NULL

#' @rdname data
#' @importFrom SingleCellExperiment cbind SingleCellExperiment
#' @importFrom stats rnbinom runif
#' @export
.mock_sc <- function(ng = 200, nc = 50, nt = 3) {
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
#' @importFrom Matrix rowSums
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
.mock_sp <- function(x, ns = 100) {
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
#' @importFrom Matrix colSums rowSums
#' @importFrom SingleCellExperiment counts
#' @importFrom stats aggregate
#' @export
.get_mgs <- function(x, n_top = 10) {
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
        .[o, ][seq_len(ifelse(n_top > nrow(.), nrow(.), n_top)), ]
    })
    z <- do.call(rbind, z)
    rownames(z) <- NULL
    return(z)
}
