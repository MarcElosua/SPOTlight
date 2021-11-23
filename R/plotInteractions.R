#' @rdname plotInteractions
#' @title Plot group interactions
#' 
#' @aliases plotHeatmap plotNetwork
#' 
#' @description This function takes in a matrix with the predicted proportions
#'   for each spot and returns a heatmap \code{which = plotHeatmap} or a network
#'    graph \code{which = plotNetwork} to show which cells are interacting
#'    spatially.
#'
#' @param x numeric matrix with rows = samples and columns = groups.
#'   Must have at least one row and column, and at least two columns.
#' @param which character string specifying the type of  
#'   visualization: one of "heatmap" or "network".
#' @param min_prop scalar specifying the value above which 
#'   a group is considered to be contributing to a given sample. 
#'   An interaction between groups i and j is counted for sample s
#'   only when both x[s, i] and x[s, j] fall above \code{min_prop}.
#' @param ... additional graphical parameters passed 
#'   to \code{plot.igraph} when \code{which = "network"}
#'   (see \code{?igraph.plotting}).
#'
#' @return base R plot
#' 
#' @author Marc Elosua Bayes & Helena L Crowell
#' 
#' @examples
#' mat <- replicate(10, rnorm(100, runif(1, -1, 1)))
#' plotNetwork(mat)
#' 
#' # specify node names
#' nms <- letters[seq_len(ncol(mat))]
#' plotNetwork(mat, vertex.label = nms)
#' 
#' # or set column names instead
#' colnames(mat) <- nms
#' plotNetwork(mat)
#' 
#' # pass additional graphical parameters for aesthetics
#' plotNetwork(mat, 
#'   edge.color = "black",
#'   vertex.color = "pink",
#'   vertex.label.font = 2,
#'   vertex.label.color = "maroon")
#'   
#' @export

plotInteractions <- function(x, 
    which = c("heatmap", "network"), 
    min_prop = 0, ...) 
{
    # check validity of input arguments
    which <- match.arg(which)   
    stopifnot(
        is.matrix(x), is.numeric(x), 
        all(dim(x) > 0), ncol(x) > 1,
        is.numeric(min_prop), length(min_prop) == 1)
        
    # get interactions table
    if (is.null(colnames(x))) 
        colnames(x) <- seq_len(ncol(x))
    df <- .count_interactions(x, min_prop)
    
    switch(which,
        heatmap = .plot_heatmap(x, df),
        network = .plot_network(x, df, ...))
}

# TODO why do we need these 2 functions?
#' @rdname plotInteractions
#' @export
plotHeatmap<- function(x, min_prop = 0, ...) 
    plotInteractions(x, "heatmap", min_prop, ...)

#' @rdname plotInteractions
#' @export
plotNetwork <- function(x, min_prop = 0, ...)
    plotInteractions(x, "network", min_prop, ...)

#' @importFrom matrixStats rowAlls
#' @importFrom utils combn
.count_interactions <- function(x, min_prop)
{
    # for each pair of groups count how many 
    # samples have value above 'min_prop'
    x <- x > min_prop
    ij <- combn(colnames(x), 2)
    y <- apply(ij, 2, \(.) sum(rowAlls(x[, ., drop = FALSE])))
    
    # construct 'data.frame'
    df <- data.frame(t(ij), y)
    names(df) <- c("from", "to", "n")
    return(df)
}

#' @import ggplot2
#' @importFrom Matrix colSums
.plot_heatmap <- function(x, df) 
{
    # assure are properly ordered
    df$from <- factor(df$from, colnames(x))
    df$to <- factor(df$to, rev(colnames(x)))
    
    # compute proportion of samples that have all groups
    i <- match(df$from, colnames(x))
    df$t <- colSums(x > 0)[i]
    df$p <- df$n / df$t
    
    ggplot(df, aes_string("from", "to", fill = "p")) +
        geom_tile() +
        coord_fixed(expand = FALSE) +
        scale_fill_viridis_c(limits = c(0, NA)) +
        theme_linedraw() +
        theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
}

#' @importFrom igraph graph_from_data_frame plot.igraph
.plot_network <- function(x, df, ...) {
    w <- scale(df$n, 1)
    g <- graph_from_data_frame(df,
        vertices = colnames(x),
        directed = FALSE)
    plot.igraph(g, edge.width = w, ...)
}
