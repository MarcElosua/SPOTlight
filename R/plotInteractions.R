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
#' @param metric character string specifying which metric to show:
#'   one of "prop" or "jaccard".
#' @param ... additional graphical parameters passed
#'   to \code{plot.igraph} when \code{which = "network"}
#'   (see \code{?igraph.plotting}).
#'
#' @return base R plot
#'
#' @author Marc Elosua Bayes & Helena L Crowell
#'
#' @examples
#' library(ggplot2)
#' mat <- replicate(8, rnorm(100, runif(1, -1, 1)))
#' # Basic example
#' plotInteractions(mat)
#'
#' ### heatmap ###
#' # This returns a ggplot object that can be modified as such
#' plotInteractions(mat, which = "heatmap") +
#'     scale_fill_gradient(low = "#f2e552", high = "#850000") +
#'     labs(title = "Interaction heatmap", fill = "proportion")
#'         
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
    metric = c("prop", "jaccard"),
    min_prop = 0, ...) {
    # check validity of input arguments
    which <- match.arg(which)
    metric <- match.arg(metric)
    stopifnot(
        is.matrix(x), is.numeric(x),
        all(dim(x) > 0), ncol(x) > 1,
        is.numeric(min_prop), length(min_prop) == 1)

    # get interactions table
    if (is.null(colnames(x))) {
        colnames(x) <- seq_len(ncol(x))
    }
    df <- .count_interactions(x, min_prop)
    df <- .statistics_interaction(x, df)

    switch(which,
        heatmap = .plot_heatmap(x, df, metric),
        network = .plot_network(x, df, metric, ...))
}

#' @importFrom matrixStats rowAlls
.count_interactions <- function(x, min_prop) {
    # for each pair of groups count how many
    # samples have value above 'min_prop'
    x <- x > min_prop
    ij <- combn(colnames(x), 2)
    y <- apply(ij, 2, function(.) sum(rowAlls(x[, ., drop = FALSE])))

    # construct 'data.frame'
    df <- data.frame(t(ij), y)
    names(df) <- c("from", "to", "n")
    
    # assure are properly ordered
    y <- colnames(x)
    df$i <- factor(df$from, y)
    df$j <- factor(df$to, rev(y))
    
    return(df)
}

.statistics_interaction <- function(x, df) {
    # compute proportion of samples that have all groups
    y <- colnames(x)
    t <- colSums(x > 0)
    i <- match(df$from, y)
    j <- match(df$to, y)
    df$ti <- t[i]
    df$tj <- t[j]
    df$pi <- df$n / df$ti
    df$pj <- df$n / df$tj
    # As suggested by @astrid12345
    # https://github.com/MarcElosua/SPOTlight/issues/42
    df$jaccard <- df$n / (df$ti + df$tj - df$n)
    return(df)
}
#' @import ggplot2
#' @importFrom Matrix colSums
.plot_heatmap <- function(x, df, metric) {
    
    # Initialize ggplot
    p <- ggplot(df)
    
    # Add pertinent layers
    if (metric == "prop") {

        # Add tile layers
        p <- p + geom_tile(aes(.data$i, .data$j, fill = .data$pi)) +
            geom_tile(aes(.data$j, .data$i, fill = .data$pj))
            
    } else if (metric == "jaccard") {
        # Add tile layers - Jaccard
        p <- p + geom_tile(aes(.data$i, .data$j, fill = .data$jaccard))
    }

    # Prettify the plot :)
    p +
        scale_fill_viridis_c("proportion", limits = c(0, NA)) +
        scale_y_discrete(limits = function(.) rev(.)) +
        coord_fixed(expand = FALSE) +
        labs(x = "From", y = "To", fill = "Proportion") +
        theme_linedraw() +
        theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
}

.plot_network <- function(x, df, metric, ...) {
    # Check necessary packages are installed and if not STOP
    .test_installed("igraph")
    
    w <- switch(metric,
        prop = scale(df[, "n"], 1),
        jaccard = df[, "jaccard"])
    
    g <- igraph::graph_from_data_frame(df,
        vertices = colnames(x),
        directed = FALSE)
    
    igraph::plot.igraph(g, edge.width = w, ...)
}
