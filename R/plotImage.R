#' @rdname plotImage
#' @name test
#' @title Plot JP(E)G/PNG/Raster/RGB images
#'
#' @description This function takes in an image-related object - path to
#'   JP(E)G/PNG file, raster object, RGBarray. It returns a ggplot object with
#'   the selected image.
#'
#' @param x A variety of objects can be passed: character string corresponding
#'   to an image file path, valid file types are JPG, JPEG and PNG. It can also
#'   take as input objects of class raster and RGB arrays. It can also take
#'   a SpatialExperiment or Seurat object from which the image will be extracted.
#' @param slice Character string indicating which image slice to use when
#'   SpatialExperiment or Seurat objects are passed. By default uses the first
#'   slice available.
#' @param alpha NOT IMPLEMENTED - single numeric between 0 and 1 determining the
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
          function(x, ..., slice = Images(x)[1]) {
            # Stop if there are no images or the name selected doesn't exist
            stopifnot(
              !is.null(Images(x)),
              slice %in% Images(x))
            x <- GetImage(x, image = slice)
            plotImage(x)
          })

#' @rdname plotImage
#' @importFrom SpatialExperiment imgRaster getImg imgData
#' @export
setMethod("plotImage", "SpatialExperiment",
          function(x, ..., slice = imgData(x)[1, "sample_id"]) {
            # Stop if there are no images or the name selected doesn't exist
            stopifnot(
              !is.null(getImg(x)),
              slice %in% imgData(spe)[1, "sample_id"])

            x <- imgRaster(spe, sample_id = slice)
            plotImage(x)
          })

#' @rdname plotImage
#' @importFrom grid unit rasterGrob
#' @export
setMethod("plotImage", "array",
          function(x) {
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
          function(x, ..., slice = NULL) {
            ggplot() +
              annotation_custom(
                grob = x,
                xmin = 0,
                xmax = ncol(x$raster),
                ymin = 0,
                ymax = nrow(x$raster)
                ) +
              coord_fixed(
                # ratio = 1
                xlim = c(0, ncol(x$raster)),
                ylim = c(0, nrow(x$raster))
                ) +
              theme_void()
          })

#' @rdname plotImage
#' @export
setMethod("plotImage", "ANY",
          function(x) {
            stop("See ?plotImage for valid image inputs")
          })
