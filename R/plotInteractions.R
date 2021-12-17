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
#' # Basic example
#' plotInteractions(mat)
#' 
#' ### heatmap ###
#' # This returns a ggplot object that can be modified as such
#' plotInteractions(mat, which = "heatmap") +
#'     scale_fill_gradient(low = "#f2e552", high = "#850000") +
#'         labs(
#'              title = "Interaction heatmap",
#'              fill = "proportion")
#' ### Network ###
#' # specify node names
#' nms <- letters[seq_len(ncol(mat))]
#' plotInteractions(mat, which = "network", vertex.label = nms)
#' 
#' # or set column names instead
#' colnames(mat) <- nms
#' plotInteractions(mat, which = "network")
#' 
#' # pass additional graphical parameters for aesthetics
#' plotInteractions(mat,
#'     which = "network",
#'     edge.color = "cyan",
#'     vertex.color = "pink",
#'     vertex.label.font = 2,
#'     vertex.label.color = "maroon")
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

#' @importFrom matrixStats rowAlls
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
    j <- match(df$to, colnames(x))
    df$t_from <- colSums(x > 0)[i]
    df$t_to <- colSums(x > 0)[j]
    df$p_from <- df$n / df$t_from
    df$p_to <- df$n / df$t_to
    
    ggplot(df, aes_string("from", "to", fill = "p")) +
        geom_tile() +
        coord_fixed(expand = FALSE) +
        scale_fill_viridis_c(limits = c(0, NA)) +
        theme_linedraw() +
        theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
}

.plot_network <- function(x, df, ...) {
    # Check necessary packages are installed and if not STOP
    .test_installed("igraph")
    
    w <- scale(df$n, 1)
    g <- igraph::graph_from_data_frame(df,
        vertices = colnames(x),
        directed = FALSE)
    igraph::plot.igraph(g, edge.width = w, ...)
}
