#' @rdname plotCorrelationMatrix
#' @name plotCorrelationMatrix
#' @title Plot Correlation Matrix
#'
#' @description This function takes in a matrix with the predicted proportions
#'   for each spot and returns a correlation matrix between cell types.
#'
#' @param x numeric matrix with rows = samples and columns = cell types
#'   Must have at least two rows and two columns.
#' @param cor.method Method to use for correlation:
#'   c("pearson", "kendall", "spearman"). By default pearson.
#' @param insig character, specialized insignificant correlation coefficients,
#'   "pch", "blank" (default). If "blank", wipe away the corresponding glyphs;
#'   if "pch", add characters (see pch for details) on corresponding glyphs.
#' @param colors character vector with three colors indicating the lower, mid,
#'   and high color. By default c("#6D9EC1", "white", "#E46726").
#' @param hc.order logical value. If TRUE, correlation matrix will be
#'   hc.ordered using hclust function.
#' @param p.mat logical value. If TRUE (default), correlation significance
#'   will be used. If FALSE arguments sig.level, insig, pch, pch.col,
#'   pch.cex are invalid.
#' @param ... additional graphical parameters passed to \code{ggcorrplot}.
#'
#' @return \code{ggplot} object
#'
#' @author Marc Elosua Bayes & Helena L Crowell
#'
#' @examples
#' set.seed(321)
#' x <- replicate(m <- 25, runif(10, 0, 1))
#' rownames(x) <- paste0("spot", seq_len(nrow(x)))
#' colnames(x) <- paste0("type", seq_len(ncol(x)))
#'
#' # The most basic example
#' plotCorrelationMatrix(x = x)
#'
#' # Showing the non-significant correlatinos
#' plotCorrelationMatrix(x = x, insig = "pch")
#'
#' # A more elaborated
#' plotCorrelationMatrix(
#'     x = x,
#'     hc.order = FALSE,
#'     type = "lower",
#'     outline.col = "lightgrey",
#'     method = "circle",
#'     colors = c("#64ccc9", "#b860bd", "#e3345d"))
#'
NULL

#' @rdname plotCorrelationMatrix
#' @importFrom Matrix colSums
#' @importFrom stats cor median
#' @import ggplot2
#' @export

plotCorrelationMatrix <- function(
    x,
    cor.method = c("pearson", "kendall", "spearman"),
    insig = c("blank", "pch"),
    colors = c("#6D9EC1", "white", "#E46726"),
    hc.order = TRUE,
    p.mat = TRUE,
    ...) {
        # Check necessary packages are installed and if not STOP
        .test_installed("ggcorrplot")
        # If the following are left undefined select
        # the first element of the vector
        cor.method <- match.arg(cor.method)
        insig <- match.arg(insig)

        stopifnot(
            is.matrix(x), is.numeric(x),
            all(dim(x) > 0), ncol(x) > 1,
            is.character(colors), length(colors) == 3,
            is.logical(hc.order), length(hc.order) == 1,
            is.logical(p.mat), length(p.mat) == 1)

        # Remove columns that are all 0
        x <- x[, colSums(x) > 0]
        corr <- cor(x)

        # Compute correlation P-value
        p.mat <- if (p.mat) {
            ggcorrplot::cor_pmat(
                x = x,
                conf_int = 0.95,
                method = cor.method)
        }

        # Plot correlation matrix as a heatmap
        ggcorrplot::ggcorrplot(
            corr = corr,
            p.mat = p.mat,
            hc.order = hc.order,
            insig = insig,
            lab = FALSE,
            colors = colors,
            ...) +
            theme(
                plot.title = element_text(hjust = 0.5, face = "bold"),
                axis.text.x = element_text(angle = 60, vjust = 1),
                axis.text = element_text(vjust = 0.5))
    }
