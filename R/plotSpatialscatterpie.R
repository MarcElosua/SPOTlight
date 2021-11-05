#' @rdname plotSpatialScatterpie
#' @name plotSpatialScatterpie
#' @title Spatial scatterpie
#'
#' @description This function takes in the coordinates of the spots and the
#'   proportions of the cell types within each spot. It returns a plot with
#'   each spot as a piechart of the proportions of the cell type composition.
#'
#' @param x Object containig the spots coordinates, it can be an object of class
#'   SpatialExperiment, Seurat, dataframe or matrix. For the latter two
#'   rownames should have the spot barcodes to match x.
#' @param y Matrix or dataframe containing the deconvoluted spots. rownames
#'   hould have the spot barcodes to match x.
#' @param img Logical TRUE or FALSE indicating whether to plot the image or not.
#'   Objects of classes accepted by \code{plotImage} can also be passed and
#'   that image will be used. By default FALSE.
#' @param slice Character string indicating which slice to plot if img is TRUE.
#'   By default uses the first image.
#' @param cell_types Vector of cell type names to plot. By default uses the
#'   column names of y.
#' @param scatterpie_alpha Numeric scalar to set the alpha of the pie charts.
#'   By default 1.
#' @param pie_scaleNumeric scalar to set the size of the pie charts.
#'   By default 0.4.
#' @return \code{ggplot} object
#'
#' @examples
#' # Filename
#' plotSpatialScatterpie("~/packages/SPOTlight/inst/extdata/SPOTlight.png")
#' path <- paste0(system.file(package="SPOTlight"), "/spotlight_ls_anterior.RDS")
#' y <- readRDS(path)[[2]]
#'
#' # array
#' # Seurat Object
#' x <- LoadData("stxBrain", type = "anterior1")
#' plotSpatialScatterpie(x, y)
#' # SpatialExperiment
#' example(read10xVisium, echo = FALSE)
#' plotSpatialScatterpie(spe, y)
NULL

#' @rdname plotSpatialScatterpie
#' @import ggplot2
#' @importFrom scatterpie geom_scatterpie
#' @export
setMethod("plotSpatialScatterpie", c("matrix", "matrix"),
          function(x, y, ...,
                   cell_types = colnames(y),
                   img = FALSE,
                   scatterpie_alpha = 1,
                   pie_scale = 0.4) {
            # Stop if x and y don't have the same number of columns or if the
            # rownames are not common between them
            # TODO add checks for image - NULL
            stopifnot(
              nrow(x) == nrow(y),
              rownames(x) %in% rownames(y)
            )
            # If image is passed add it as the base layer, if not, no image
            # Need to use isFALSE bc img can have many different inputs
            # Set ymax to overlap image and piecharts
            if (isFALSE(img)) {
              p <- ggplot()
              ymax <- 0
            } else {
              p <- plotImage(x = img) + scale_y_reverse()
              ymax <- max(p$coordinates$limits$y)
            }

            # TODO Step below, in the new implementation instead of cbind use merge
            # We need to make sure the obtained matrices have col and rownames
            # df <- data.frame(merge(x, y, by = "row.names"))
            df <- data.frame(cbind(x, y))

            # Plot
            p + geom_scatterpie(
              data = df,
              aes(x = imagecol,
                  y = abs(imagerow - ymax)),
              cols = cell_types,
              color = NA,
              alpha = scatterpie_alpha,
              pie_scale = pie_scale) +
              coord_fixed() +
              theme_void()
          })

#' @rdname plotSpatialScatterpie
#' @importFrom Seurat GetTissueCoordinates Images GetImage
#' @export
setMethod("plotSpatialScatterpie", c("Seurat", "ANY"),
          function(x, y, ..., slice = Images(x)[1], img = FALSE) {
            # Stop if there image is to be extracted from Seurat object but
            # no images present or slice selected doesn't exist
            if (isTRUE(img)) {
              stopifnot(
                !is.null(Images(x)),
                slice %in% Images(x))
            }

            # If img is NULL extract image from Seurat object
            if (img) img <- GetImage(x, image = slice)

            ## Extract spatial coordinates
            x <- as.matrix(GetTissueCoordinates(x, image = slice))

            plotSpatialScatterpie(x, y, img = img)
          })

#' @rdname plotSpatialScatterpie
#' @importFrom SpatialExperiment spatialCoords imgData
#' @export
setMethod("plotSpatialScatterpie", c("SpatialExperiment", "ANY"),
          function(x, y, ..., slice = imgData(x)[1, "sample_id"], img = FALSE) {

            # TODO Stop if there are no images or the name selected doesn't exist
            stopifnot(
              isTRUE(img) & !is.null(getImg(x)),
              isTRUE(img) & slice %in% imgData(spe)[1, "sample_id"])

            # If img is NULL extract image from Seurat object
            if (img) img <- imgRaster(spe, sample_id = slice)

            ## Extract spot barcodes
            barcodes <- colnames(x)
            ## Extract spatial coordinates
            x <- as.matrix(spatialCoords(x)[, c("x", "y")])
            ## Add barcodes to coord matrix
            rownames(x) <- barcodes
            colnames(x) <- c("imagecol", "imagerow")

            plotSpatialScatterpie(x, y)
          })

#' @rdname plotSpatialScatterpie
#' @export
setMethod("plotSpatialScatterpie", c("data.frame", "ANY"),
          function(x, y, ...) {

            # Convert data.frame to matrix
            x <- as.matrix(x)

            plotSpatialScatterpie(x, y)
          })

#' @rdname plotSpatialScatterpie
#' @export
setMethod("plotSpatialScatterpie", c("ANY", "data.frame"),
          function(x, y, ...) {
            # Convert data.frame to matrix
            y <- as.matrix(y)

            plotSpatialScatterpie(x, y)
          })

#' @rdname plotSpatialScatterpie
#' @export
setMethod("plotSpatialScatterpie", c("ANY", "ANY"),
          function(x) {
            stop("See ?plotSpatialScatterpie for valid inputs")
          })
