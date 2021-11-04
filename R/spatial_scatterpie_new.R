#' @rdname plotSpatialscatterpie
#' @name test
#' @title Spatial scatterpie
#'
#' @description ...
#'
#' @param x Object containig the spots coordinates
#' @param y Matrix containing the deconvoluted spots
#'
#' @return \code{ggplot} object
#'
#' @examples
#' # Filename
#' plotSpatialscatterpie("~/packages/SPOTlight/inst/extdata/SPOTlight.png")
#' path <- paste0(system.file(package="SPOTlight"), "/spotlight_ls_anterior.RDS")
#' y <- readRDS(path)[[2]]
#' plotSpatialscatterpie(path)
#' # array
#' # Seurat Object
#' # SpatialExperiment
NULL

img <- so
x <- spatial_coord[, c("imagerow", "imagecol")]
#' @rdname plotSpatialscatterpie
#' @import ggplot2
#' @importFrom scatterpie geom_scatterpie
#' @export
setMethod("plotSpatialscatterpie", "matrix", "matrix",
          function(x, y, ..., img = NULL) {

            # If image is passed add it as the base layer, if not, leave it empty
            if (!is.null(img)) {
              p <- plotImage(x = img)
            } else p <- ggplot()

            df <- cbind(x, y)

            # Get ymax to overlap image and piecharts
            ymax <- ifelse(is.null(img), 0, max(p$coordinates$limits$y))
            p + geom_scatterpie(
                data = df,
                aes(x = imagecol,
                    y = abs(imagerow - ymax)),
                cols = cell_types_all,
                color = NA,
                alpha = scatterpie_alpha,
                pie_scale = pie_scale) +
                { if (is.null(img)) scale_y_reverse() } +
                { if (is.null(img)) coord_fixed(1) } +
                theme_void()

          })

#' @rdname plotSpatialscatterpie
#' @importFrom Seurat GetTissueCoordinates Images
#' @export
setMethod("plotSpatialscatterpie", "Seurat", "ANY",
          function(x, y, ..., slice = Images(x)[1]) {
            # Stop if there are no images or the name selected doesn't exist
            stopifnot(
              !is.null(Images(x)),
              slice %in% Images(x))

            ## Extract spatial coordinates
            x <- as.matrix(GetTissueCoordinates(x, image = slice))

            plotSpatialscatterpie(x, y)

          })

#' @rdname plotSpatialscatterpie
#' @importFrom SpatialExperiment spatialCoords imgData
#' @export
setMethod("plotSpatialscatterpie", "SpatialExperiment", "ANY",
          function(x, y, ..., slice = imgData(x)[1, "sample_id"]) {
            ## Extract spot barcodes
            barcodes <- colnames(x)
            ## Extract spatial coordinates
            x <- as.matrix(spatialCoords(x)[, c("x", "y")])
            ## Add barcodes to coord matrix
            rownames(x) <- barcodes
            colnames(x) <- c("imagecol", "imagerow")
            plotSpatialscatterpie(x, y)
          })

#' @rdname plotSpatialscatterpie
#' @importFrom SpatialExperiment spatialCoords imgData
#' @export
setMethod("plotSpatialscatterpie", "data.frame", "ANY",
          function(x, y, ...) {
            # Convert data.frame to matrix
            x <- as.matrix(x)

            plotSpatialscatterpie(x, y)
          })

#' @rdname plotSpatialscatterpie
#' @export
setMethod("plotSpatialscatterpie", "ANY", "ANY",
          function(x) {
            stop("See ?plotSpatialscatterpie for valid image inputs")
          })
