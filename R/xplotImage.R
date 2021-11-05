#' @rdname plotImage
#' @name test
#' @title Plot JP(E)G/PNG image
#' 
#' @description ...
#' 
#' @param x character string corresponding to an image file path.
#'   Valid file types are JPG, JPEG and PNG.
#' @param alpha single numeric in between 0 and 1 determining the
#'   image opacity. Lower values correspond to more transparency.
#'
#' @return \code{ggplot} object
#'
#' @examples
#' plotImage("~/packages/SPOTlight/inst/extdata/SPOTlight.png")
NULL

#' @rdname plotImage
#' @importFrom jpeg readJPEG
#' @importFrom png readPNG
#' @export
setMethod("plotImage", "character", 
    function(x) {
        stopifnot(file.exists(x))
        typ <- c("jpg", "jpeg", "png")
        pat <- paste0(".", typ, "$")
        idx <- vapply(pat, grepl, x = x, logical(1))
        if (!any(idx)) stop("'x' should be JPG, JPEG or PNG")
        x <- switch(typ[idx], png = readPNG(x), readJPEG(x))
        plotImage(x)
    })

#' @rdname plotImage
#' @import ggplot2
#' @importFrom grid unit rasterGrob
#' @export
setMethod("plotImage", "array", 
    function(x) {
        x <- rasterGrob(x, 
            interpolate = FALSE, 
            width = unit(1, "npc"),
            height = unit(1, "npc"))
        ggplot() + 
            annotation_custom(x) +
            coord_fixed(
                xlim = c(0, ncol(x$raster)), 
                ylim = c(0, nrow(x$raster))) + 
            theme_void()
    })
