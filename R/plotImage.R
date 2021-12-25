#' @rdname plotImage
#' @name plotImage
#' @title Plot JP(E)G/PNG/Raster/RGB images
#'
#' @description This function takes in an image-related object - path to
#'   JP(E)G/PNG file, raster object, RGBarray. It returns a ggplot object with
#'   the selected image.
#'
#' @param x A variety of objects can be passed: character string corresponding
#'   to an image file path, valid file types are JPG, JPEG and PNG. It can also
#'   take as input objects of class raster and RGB arrays. It can also take
#'   a SpatialExperiment or Seurat object from which the
#'   image will be extracted.
#' @param slice Character string indicating which image slice to use when
#'   SpatialExperiment or Seurat objects are passed. By default uses the first
#'   slice available.
#' @param alpha NOT IMPLEMENTED - single numeric between 0 and 1 determining the
#'   image opacity. Lower values correspond to more transparency.
#' @return \code{ggplot} object
#'
#' @author Marc Elosua Bayes & Helena L Crowell
#'
#' @examples
#' # Filename
#' path <- paste0(system.file(package = "SPOTlight"), "/extdata/SPOTlight.png")
#' plotImage(x = path)
#' # array
#' png_img <- png::readPNG(path)
#' plotImage(png_img)
#' # Seurat Object
#' # library(SeuratData)
#' # so <- LoadData("stxBrain", type = "anterior1")
#' # plotImage(so)
#' # SpatialExperiment
#' library(TENxVisiumData)
#' spe <- MouseKidneyCoronal()
#' plotImage(spe)
NULL

#' @rdname plotImage
#' @export
setMethod(
    "plotImage", "character",
    function(x) {
        # Check necessary packages are installed and if not STOP
        .test_installed(c("jpeg", "png"))

        stopifnot(file.exists(x))
        typ <- c("jpg", "jpeg", "png")
        pat <- paste0(".", typ, "$")
        idx <- vapply(pat, grepl, x = x, logical(1))
        if (!any(idx)) {
            stop("'x' should be of file type JPG, JPEG or PNG")
        }
        x <- switch(typ[idx],
            png = png::readPNG(x),
            jpeg::readJPEG(x)
        )
        plotImage(x)
    }
)

#' @rdname plotImage
#' @importFrom SeuratObject GetImage Images
#' @export
setMethod(
    "plotImage", "Seurat",
    function(x, slice = Images(x)[1]) {
        # Stop if there are no images or the name selected doesn't exist
        stopifnot(
            !is.null(Images(x)),
            slice %in% Images(x)
        )
        # Extract Image in raster format
        x <- GetImage(x, image = slice, mode = "raster")
        # pass as a matrix
        plotImage(as.matrix(x))
    }
)

#' @rdname plotImage
#' @export
setMethod(
    "plotImage", "SpatialExperiment",
    function(x, slice = imgData(x)[1, "sample_id"]) {
        .test_installed(c("SpatialExperiment"))

        # Stop if there are no images or the name selected doesn't exist
        stopifnot(
            !is.null(SpatialExperiment::getImg(x)),
            slice %in% SpatialExperiment::imgData(x)[1, "sample_id"]
        )
        # Convert to raster
        x <- SpatialExperiment::imgRaster(x, sample_id = slice)
        # pass as a matrix
        plotImage(as.matrix(x))
    }
)

#' @rdname plotImage
#' @import ggplot2
#' @export
setMethod(
    "plotImage", "array",
    function(x, alpha = NULL) {
        # Check necessary packages are installed and if not STOP
        .test_installed("grid")

        x <- grid::rasterGrob(x,
            interpolate = FALSE,
            width = grid::unit(1, "npc"),
            height = grid::unit(1, "npc")
        )

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
    }
)

#' @rdname plotImage
#' @export
setMethod(
    "plotImage", "ANY",
    function(x) stop("See ?plotImage for valid image inputs")
)
