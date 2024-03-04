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
#'   a SpatialExperiment from which the image will be extracted.
#' @param slice Character string indicating which image slice to use when
#'   SpatialExperiment objects are passed. By default uses the first
#'   slice available.
#' @return \code{ggplot} object
#'
#' @author Marc Elosua Bayes & Helena L Crowell
#'
#' @examples
#' # Filename
#' path <- file.path(
#'   system.file(package = "SPOTlight"), 
#'   "extdata/SPOTlight.png")
#' plotImage(x = path)
#' # array
#' png_img <- png::readPNG(path)
#' plotImage(png_img)
#' # SpatialExperiment
NULL
#' @export
plotImage <- function(x, slice = NULL) {
    # check validity of input arguments
    stopifnot(
        # Check for valid x classes
        is.matrix(x) | is.character(x) | is.array(x) | is(x, "rastergrob") |
            is(x, "SpatialExperiment"),
        # Check for valid slice classes
        is.null(slice) | is.character(slice))
    
    if (!is.array(x))
        x <- .extract_image(x)
        
    # Plot image
    plt <- .plot_image(x)
}
