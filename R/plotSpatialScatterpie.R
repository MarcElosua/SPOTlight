#' @rdname plotSpatialScatterpie
#' @name plotSpatialScatterpie
#' @title Spatial scatterpie
#'
#' @description This function takes in the coordinates of the spots and the
#'   proportions of the cell types within each spot. It returns a plot where
#'   each spot is a piechart showing proportions of the cell type composition.
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
#' @author Marc Elosua Bayes & Helena L Crowell
#'
#' @examples
#' set.seed(321)
#' # Coordinates
#' x <- replicate(2, rnorm(100))
#' rownames(x) <- paste0("spot", seq_len(nrow(x)))
#' colnames(x) <- c("imagecol", "imagerow")
#' # Proportions
#' y <- replicate(m <- 5, runif(nrow(x), 0, 1))
#' y <- prop.table(y, 1)
#' rownames(y) <- paste0("spot", seq_len(nrow(y)))
#' colnames(y) <- paste0("type", seq_len(ncol(y)))
#' (plt <- plotSpatialScatterpie(x = x, y = y))
NULL
#' @rdname plotSpatialScatterpie
#' @import ggplot2
#' @export
setMethod(
    "plotSpatialScatterpie",
    c("matrix", "matrix"),
    function(x, y, ...,
    cell_types = colnames(y),
    img = FALSE,
    scatterpie_alpha = 1,
    pie_scale = 0.4) {
        # Check necessary packages are installed and if not STOP
        .test_installed("scatterpie")

        # Stop if x and y don't have the same number of columns or if the
        # rownames are not common between them
        # TODO add checks for image - NULL
        stopifnot(
            nrow(x) == nrow(y),
            all(rownames(x) %in% rownames(y))
        )

        # If image is passed add it as the base layer, if not, no image
        # Need to use isFALSE bc img can have many different inputs
        # Set ymax to overlap image and piecharts
        if (isFALSE(img)) {
            p <- ggplot() +
                coord_fixed()
            ymax <- 0
        } else {
            p <- plotImage(x = img) + scale_y_reverse()
            ymax <- max(p$coordinates$limits$y)
        }

        # merge by row names (by=0 or by="row.names")
        df <- merge(x, y, by = 0, all = TRUE)

        # Plot
        p + scatterpie::geom_scatterpie(
            data = df,
            aes(
                x = imagecol,
                y = abs(imagerow - ymax)
            ),
            cols = cell_types,
            color = NA,
            alpha = scatterpie_alpha,
            pie_scale = pie_scale
        ) +
            # Below not needed bc comes from plotImage
            # coord_fixed() +
            theme_void()
    }
)

#' @rdname plotSpatialScatterpie
#' @importFrom SeuratObject GetTissueCoordinates GetImage Images
#' @export
setMethod(
    "plotSpatialScatterpie",
    c("Seurat", "ANY"),
    function(x, y, ...,
    slice = Images(x)[1],
    img = FALSE) {
        # Stop if there image is to be extracted from Seurat object but
        # no images present or slice selected doesn't exist
        if (isTRUE(img)) {
              stopifnot(
                  !is.null(Images(x)),
                  slice %in% Images(x)
              )
          }

        # If 'img = TRUE' extract image from Seurat object
        if (img) img <- GetImage(x, image = slice)

        # Extract spatial coordinates
        x <- as.matrix(GetTissueCoordinates(x, image = slice))

        plotSpatialScatterpie(x, y, img = img)
    }
)

#' @rdname plotSpatialScatterpie
#' @export
setMethod(
    "plotSpatialScatterpie",
    c("SpatialExperiment", "ANY"),
    function(x, y, ...,
    slice = imgData(x)[1, "sample_id"],
    img = FALSE) {
        # Check necessary packages are installed and if not STOP
        .test_installed("SpatialExperiment")

        # TODO Stop if there are no images or the name selected doesn't exist
        if (isTRUE(img)) {
              stopifnot(
                  !is.null(SpatialExperiment::getImg(x)),
                  slice %in% SpatialExperiment::imgData(spe)[1, "sample_id"]
              )
          }

        # If 'img = TRUE' extract image from Seurat object
        if (img) img <- SpatialExperiment::imgRaster(spe, sample_id = slice)

        ## Extract spot barcodes
        barcodes <- colnames(x)

        ## Extract spatial coordinates
        x <- as.matrix(SpatialExperiment::spatialCoords(x)[, 1:2])

        ## Add barcodes to coord matrix & change colnames
        rownames(x) <- barcodes
        colnames(x) <- c("imagecol", "imagerow")

        plotSpatialScatterpie(x, y)
    }
)

#' @rdname plotSpatialScatterpie
#' @export
setMethod(
    "plotSpatialScatterpie", c("data.frame", "ANY"),
    function(x, y, ...) plotSpatialScatterpie(as.matrix(x), y)
)

#' @rdname plotSpatialScatterpie
#' @export
setMethod(
    "plotSpatialScatterpie", c("ANY", "data.frame"),
    function(x, y, ...) plotSpatialScatterpie(x, as.matrix(y))
)

#' @rdname plotSpatialScatterpie
#' @export
setMethod(
    "plotSpatialScatterpie", c("ANY", "ANY"),
    function(x, y, ...) stop("See ?plotSpatialScatterpie for valid inputs")
)
