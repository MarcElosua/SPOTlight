#' @title Plot NMF topic profiles
#' @description ...
#'
#' @param x \code{\link{NMFfit}} object
#' @param y vector of group labels. Should be of length \code{ncol(coef(x))}.
#' @param facet logical indicating whether to stratify by group. 
#'   If \code{FALSE}, weights will be averaged across cells for each group.
#' @param min_prop scalar in [0,1]. When \code{facet = TRUE}, 
#'   only cells with a weight > \code{min_prop} will be included.
#' @param ncol integer scalar specifying the number of facet columns.
#'
#' @return \code{ggplot} object
#' 
#' @examples
#' x <- .mock_sc()
#' y <- .mock_sp(x)
#' z <- .get_mgs(x)
#' res <- SPOTlight(x, y, groups = x$type, mgs = z, group_id = "type")
#' 
#' @importFrom methods is
#' @importFrom NMF coef
#' @import ggplot2
#' 
#' @export

setMethod(
    "plotTopicProfiles", 
    c("NMF", "character"),
    function(x, y, 
        facet = FALSE, 
        min_prop = 0.1, 
        ncol = NULL)
    {
        # check validity of input arguments
        stopifnot(
            is(x, "NMF"), length(y) == ncol(coef(x)),
            is.logical(facet), length(facet) == 1,
            is.numeric(min_prop), length(min_prop) == 1, 
            min_prop >= 0, min_prop <= 1)
        
        # get group proportions
        mat <- prop.table(t(coef(x)), 1)
        
        if (TRUE) {
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
            f <- facet_wrap(~ group, ncol = "a", scales = "free_x")
        } else {
            # get topic medians
            df <- aggregate(mat, list(y), median)[, -1]
            
            # stretch for plotting
            df <- data.frame(
                weight = unlist(df),
                group = rep(colnames(df), each = nrow(df)),
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
            guides(col = guide_legend()) +
            scale_size_continuous(range = c(0, 3)) +
            scale_color_continuous(low = "lightgrey", high = "blue") +
            xlab(ifelse(facet), "", x) + 
            theme_bw() + theme(
                panel.grid = element_blank(),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1))
    })
