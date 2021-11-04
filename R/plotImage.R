#' @rdname plotImage
#' @name test
#' @title Plot JP(E)G/PNG image
#'
#' @description ...
#'
#' @param x character string corresponding to an image file path.
#'   Valid file types are JPG, JPEG and PNG.
#' @param alpha single numeric in between 0 and 1 determining the
#'   image opacity. Lower values correspond to more transparency. \TODO
#'
#' @return \code{ggplot} object
#'
#' @examples
#' # Filename
#' plotImage("~/packages/SPOTlight/inst/extdata/SPOTlight.png")
#' path <- paste0(system.file(package="SPOTlight"), "/spatial/tissue_lowres_image.png")
#' plotImage(path)
#' # array
#' png_img <- png::readPNG(path)
#' plotImage(png_img)
#' # Seurat Object
#' so <- LoadData("stxBrain", type = "anterior1")
#' plotImage(so)
#' # SpatialExperiment
#' spe0 <- as.SingleCellExperiment(LoadData("stxBrain", type = "anterior1"))
#' plotImage(spe0)
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
            # If typ[idx] returns png readPNG if not readJPEG
            x <- switch(typ[idx], png = readPNG(x), readJPEG(x))
            plotImage(x)
          })

#' @rdname plotImage
#' @importFrom Seurat GetImage Images
#' @export
setMethod("plotImage", "Seurat",
          function(x) {
            stopifnot(is.null(Images(x)))
            print("Seurat")
            x <- GetImage(x)
            plotImage(x)
          })

#' @rdname plotImage
#' @importFrom SpatialExperiment imgRaster
#' @export
setMethod("plotImage", "SpatialExperiment",
          function(x) {
            print("SpatialExperiment")
            x <- imgRaster(spe)
            plotImage(x)
          })

#' @rdname plotImage
#' @importFrom grid unit rasterGrob
#' @export
setMethod("plotImage", "array",
          function(x) {
            print("array")
            x <- rasterGrob(x,
                            interpolate = FALSE,
                            width = unit(1, "npc"),
                            height = unit(1, "npc"))
            plotImage(x)
          })

#' @rdname plotImage
#' @import ggplot2
#' @export
setMethod("plotImage", "rastergrob",
          function(x) {
            print("rastergrob")
            ggplot() +
              annotation_custom(x) +
              coord_fixed(
                xlim = c(0, ncol(x$raster)),
                ylim = c(0, nrow(x$raster))) +
              theme_void()
          })

plotImage(paste0(system.file(package = "SPOTlight"), "/allen_cortex_dwn.rds"))
