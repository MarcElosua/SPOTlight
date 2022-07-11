#' @rdname plotTopicProfiles
#' @name plotTopicProfiles
#' @title Plot NMF topic profiles
#'
#' @description This function takes in the fitted NMF model and returns the
#'   topic profiles learned for each cell \code{facet = FALSE} or cell type
#'   \code{facet = TRUE}. Ideal training will return all the cell from the same
#'   cell type to share a unique topic profile.
#'
#' @param x \code{\link{NMFfit}} object
#' @param y vector of group labels. Should be of length \code{ncol(coef(x))}.
#' @param facet logical indicating whether to stratify by group.
#'   If \code{FALSE} (default), weights will be the median across cells
#'   for each group (point = topic weight for a given cell type).
#'   If \code{TRUE}, cell-specific weights will be shown
#'   (point = topic weight of a given cell).
#' @param min_prop scalar in [0,1]. When \code{facet = TRUE},
#'   only cells with a weight > \code{min_prop} will be included.
#' @param ncol integer scalar specifying the number of facet columns.
#' 
#' @return \code{ggplot} object
#'
#' @author Marc Elosua Bayes & Helena L Crowell
#'
#' @examples
#' library(ggplot2)
#' x <- mockSC()
#' y <- mockSP(x)
#' z <- getMGS(x)
#' 
#' res <- SPOTlight(x, y,
#'     groups = x$type,
#'     mgs = z,
#'     group_id = "type",
#'     verbose = FALSE)
#'
#' plotTopicProfiles(res[[3]], x$type, facet = TRUE)
#' plotTopicProfiles(res[[3]], x$type, facet = FALSE)
NULL

#' @rdname plotTopicProfiles
#' @importFrom stats aggregate median
#' @import ggplot2
#' @export
plotTopicProfiles <- function(
    x,
    y,
    facet = FALSE,
    min_prop = 0.01,
    ncol = NULL) {
    # Convert y to character
    y <- as.character(y)
    
    # check validity of input arguments
    stopifnot(
        is(x, "list"),
        all(sort(names(x)) == sort(c("w", "d", "h", "misc"))),
        is.character(y),
        length(y) == ncol(x$h),
        setequal(
            colnames(x$w), paste0("topic_", seq_len(length(unique(y))))
            ),
        is.logical(facet), length(facet) == 1,
        is.numeric(min_prop), length(min_prop) == 1,
        is.null(ncol) | (is.numeric(ncol) & length(ncol) == 1))
    
    # get proportion of topic contribution by cell
    mat <- prop.table(t(x$h), 1)
    
    if (facet) {
        # stretch for plotting
        df <- data.frame(
            id = seq_len(nrow(mat)),
            weight = c(mat),
            group = rep(y, ncol(mat)),
            topic = rep(seq_len(ncol(mat)), each = nrow(mat)))

        # drop cells with 'weight < min_prop'
        df <- df[df$weight >= min_prop, ]

        # set aesthetics
        x <- "id"
        f <- facet_wrap(~group, ncol = ncol, scales = "free_x")
    } else {
        # get topic medians
        df <- aggregate(mat, list(y), median)[, -1]
        rownames(df) <- unique(y)
        
        # stretch for plotting
        df <- data.frame(
            weight = unlist(df),
            group = rep(rownames(df), each = nrow(df)),
            topic = rep(seq_len(nrow(df)), ncol(df)))

        # set aesthetics
        x <- "group"
        f <- NULL
    }
    # fix topic order
    df$topic <- factor(df$topic, seq_along(unique(y)))

    # render plot
    ggplot(df, aes_string(x, "topic",
        col = "weight", size = "weight")) +
        f + geom_point() +
        guides(col = guide_legend(override.aes = list(size = 2))) +
        scale_size_continuous(range = c(0, 3)) +
        scale_color_continuous(low = "lightgrey", high = "#3d2bff") +
        xlab(if (facet) x) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            legend.key.size = unit(0.5, "lines"),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1))
}

