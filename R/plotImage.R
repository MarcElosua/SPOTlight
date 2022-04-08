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
#' # Seurat Object
#' # library(SeuratData)
#' # so <- LoadData("stxBrain", type = "anterior1")
#' # plotImage(so)
#' # SpatialExperiment
NULL
#' @export
plotImage <- function(x, slice = NULL) {
    # check validity of input arguments
    stopifnot(
        # Check for valid x classes
        is.matrix(x) | is.character(x) | is.array(x) | is(x, "rastergrob") | 
            is(x, "Seurat") | is(x, "SpatialExperiment"),
        # Check for valid slice classes
        is.null(slice) | is.character(slice))
    
    
    if (is.array(x)) {
        # Plot image
        plt <- .plot_image(x)
    } else {
        # Extract image
        x <- .extract_image(x)
        # Plot image
        plt <- .plot_image(x)
    }
}

# Take an array representing an image and plot it with ggplot2
#' @import ggplot2
.plot_image <- function(x) {
    # Check necessary packages are installed and if not STOP
    .test_installed(c("grid", "ggplot2"))
    
    x <- grid::rasterGrob(x,
        interpolate = FALSE,
        width = grid::unit(1, "npc"),
        height = grid::unit(1, "npc"))
    
    ggplot() +
        annotation_custom(
            grob = x,
            xmin = 0,
            xmax = ncol(x$raster),
            ymin = 0,
            ymax = nrow(x$raster)) + 
        coord_fixed(
            xlim = c(0, ncol(x$raster)),
            ylim = c(0, nrow(x$raster))) + 
        theme_void()
}

# Extract image and convert it to array from allowed classes
.extract_image <- function(x, slice = NULL) {
    # Iterate over all the accepted classes and conver the image to array
    if (is.character(x)) {
        .test_installed(c("jpeg", "png"))
        
        # Check if the file exists
        stopifnot(file.exists(x))
        
        # Check the file is in the right format
        typ <- c("jpg", "jpeg", "png")
        pat <- paste0(".", typ, "$")
        idx <- vapply(pat, grepl, x = x, logical(1), ignore.case = TRUE)
        if (!any(idx)) {
            stop("'x' should be of file type JPG, JPEG or PNG")
        }
        
        # Read file
        x <- switch(typ[idx],
            png = png::readPNG(x),
            jpeg::readJPEG(x))
        
    } else if (is(x, "Seurat")) {
        .test_installed(c("SeuratObject"))
        # Stop if there are no images or the name selected doesn't exist
        stopifnot(
            !is.null(SeuratObject::Images(x)),
            slice %in% SeuratObject::Images(x))
        
        # If image is null use the first slice
        if (is.null(slice)) 
            slice <- SeuratObject::Images(x)[1]
        
        # Extract Image in raster format
        x <- SeuratObject::GetImage(x, image = slice, mode = "raster")
        # Conver to matrix
        x <- as.matrix(x)
        
    } else if (is(x, "SpatialExperiment")) {
        
        .test_installed(c("SpatialExperiment"))
        
        # Stop if there are no images or the name selected doesn't exist
        stopifnot(
            !is.null(SpatialExperiment::getImg(x)),
            slice %in% SpatialExperiment::imgData(x)[1, "sample_id"]
        )
        
        # If image is null use the first slice
        if (is.null(slice)) 
            slice <- SpatialExperiment::imgData(x)[1, "sample_id"]
        
        # Convert to raster
        x <- SpatialExperiment::imgRaster(x, sample_id = slice)
        x <- as.matrix(x)
    } else {
        stop("Couldn't extract image, See ?plotImage for valid image inputs.")
    }
    return(x)
}
