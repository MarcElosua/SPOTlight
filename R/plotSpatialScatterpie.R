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
#'   rownames should have the spot barcodes to match x. If a matrix it has to
#'   of dimensions nrow(y) x 2 where the columns are the x and y coordinates
#'   in that order.
#' @param y Matrix or dataframe containing the deconvoluted spots. rownames
#'   need to be the spot barcodes to match to x.
#' @param img Logical TRUE or FALSE indicating whether to plot the image or not.
#'   Objects of classes accepted by \code{plotImage} can also be passed and
#'   that image will be used. By default FALSE.
#' @param slice Character string indicating which slice to plot if img is TRUE.
#'   By default uses the first image.
#' @param cell_types Vector of cell type names to plot. By default uses the
#'   column names of y.
#' @param scatterpie_alpha Numeric scalar to set the alpha of the pie charts.
#'   By default 1.
#' @param pie_scale Numeric scalar to set the size of the pie charts.
#'   By default 0.4.
#' @param degrees From SpatialExperiment rotateImg. For clockwise (degrees > 0)
#'  and counter-clockwise (degrees < 0) rotation. By default NULL.
#' @param axis From SpatialExperiment mirrorImg. When a SpatialExperiment object
#'   is passed as the image return the mirror image. For horizontal (axis = "h")
#'    and vertical (axis = "v") mirroring. By default NULL.
#' @param ... additional parameters to geom_scatterpie
#' @return \code{ggplot} object
#'
#' @author Marc Elosua Bayes & Helena L Crowell
#'
#' @examples
#' set.seed(321)
#' 
#' # Coordinates
#' x <- replicate(2, rnorm(100))
#' rownames(x) <- paste0("spot", seq_len(nrow(x)))
#' colnames(x) <- c("imagecol", "imagerow")
#' 
#' # Proportions
#' y <- replicate(m <- 5, runif(nrow(x), 0, 1))
#' y <- prop.table(y, 1)
#' 
#' rownames(y) <- paste0("spot", seq_len(nrow(y)))
#' colnames(y) <- paste0("type", seq_len(ncol(y)))
#' 
#' (plt <- plotSpatialScatterpie(x = x, y = y))
NULL

#' @rdname plotSpatialScatterpie
#' @import ggplot2
#' @export
plotSpatialScatterpie <- function(
    x,
    y,
    cell_types = colnames(y),
    img = FALSE,
    slice = NULL,
    scatterpie_alpha = 1,
    pie_scale = 0.4,
    degrees = NULL,
    axis = NULL,
    ...) {
    # Check necessary packages are installed and if not STOP
    .test_installed("scatterpie")
    
    # Class checks
    stopifnot(
        # Check x inputs
        is.matrix(x) | is.data.frame(x) |
            is(x, "Seurat") | is(x, "SpatialExperiment"),
        # Check y inputs
        is.matrix(y) | is.data.frame(y),
        # cell_types needs to be a character with max length = ncol(y)
        is.character(cell_types) & length(cell_types) <= ncol(y),
        # Check img
        # img not checked since its checked in plotImage()
        # Check slice name
        is.character(slice) | is.null(slice),
        # Check plotting parameters are numeric
        is.numeric(scatterpie_alpha),
        is.numeric(pie_scale),
        is.numeric(degrees) | is.null(degrees),
        axis %in% c("h", "v") | is.null(axis)
    )
    
    # If image is passed add it as the base layer, if not, no image
    # Need to use isFALSE bc img can have many different inputs
    # Set ymax to overlap image and piecharts
    if (isFALSE(img)) {
        p <- ggplot() +
            coord_fixed()
        ymax <- 0
    } else {
        # Extract image from Seurat or SE objects when img is TRUE
        # If image is not TRUE and not FALSE an acceptable class for plotImage
        # has been passed
        if (is(x, "Seurat") | is(x, "SpatialExperiment") & isTRUE(img)) {
            img <- .extract_image(x, slice)
            
            # Rotate or mirror image if dots don't overlay properly
            if (is(x, "SpatialExperiment")) {
                .test_installed("SpatialExperiment")
                
                ## Rotate image if needed
                if (!is.null(degrees)) {
                    img <- SpatialImage(as.raster(img))
                    img <- as(img, "LoadedSpatialImage")
                    img <- SpatialExperiment::rotateImg(img, degrees = degrees)
                    img <- as.raster(img)
                }
                
                ## Make mirror image if necessary
                if (!is.null(axis)) {
                    img <- SpatialImage(as.raster(img))
                    img <- as(img, "LoadedSpatialImage")
                    img <- SpatialExperiment::mirrorImg(img, axis = axis)
                    img <- as.raster(img)
                }
            }
        }
        
        p <- plotImage(x = img)
        ymax <- max(p$coordinates$limits$y)
    }
    
    # Extract coordinate matrix from x
    if (!is.matrix(x))
        x <- .extract_coord(x = x, slice = slice, img = img)
    # Check colnames
    x <- .x_cnames(x)
    
    # Convert y to matrix format
    if (!is.matrix(x)) {
        y <- as.matrix(x)
    }
    
    # Stop if x and y don't have the same number of columns or if the
    # rownames are not common between them
    stopifnot(
        nrow(x) == nrow(y),
        all(rownames(x) %in% rownames(y)))
    
    # merge by row names (by=0 or by="row.names")
    df <- merge(x, y, by = 0, all = TRUE)

    # Plot
    p + scatterpie::geom_scatterpie(
        data = df,
        aes(
            x = coord_x,
            y = abs(coord_y - ymax)
        ),
        cols = cell_types,
        color = NA,
        alpha = scatterpie_alpha,
        pie_scale = pie_scale,
        ...) +
        # Below not needed bc comes from plotImage
        # coord_fixed() +
        theme_void() + 
        theme(legend.key.size = unit(0.5, "lines"))
    }

.x_cnames <- function(x) {
    # If the column names of x aren't right fix them
    cnames <- c("coord_x", "coord_y")
    if (!all(colnames(x) %in% cnames)) {
        colnames(x) <- cnames
    }
    x
}

# Coordinates and return a matrix object where each row is a spot and the
# columns are the x and y coordinates
.extract_coord <- function(x, slice, img) {
    # Iterate over all the accepted classes and return spot coordinates
    if (is.data.frame(x)) {
        # Convert to matrix
        x <- as.matrix(x)
    } else if (is(x, "Seurat")) {
        .test_installed(c("SeuratObject"))
        # Stop if there are no images or the name selected doesn't exist
        stopifnot(
            # Stop if there are no images
            !is.null(SeuratObject::Images(x)),
            # Stop if the image doesn't exist
            slice %in% SeuratObject::Images(x))
        
        # If image is null use the first slice
        if (is.null(slice) & img) 
            slice <- SeuratObject::Images(x)[1]
        
        # Extract Image in raster format
        # If 'img = TRUE' extract image from Seurat object
        # TODO Check if we can delete this since we can pass Seurat to plotImage
        # if (img) img <- GetImage(x, image = slice, mode = "raster")
        
        # Extract spatial coordinates
        x <- as.matrix(SeuratObject::GetTissueCoordinates(x, image = slice))
    } else if (is(x, "SpatialExperiment")) {
        
        .test_installed(c("SpatialExperiment"))
        
        # Stop if there are no images or the name selected doesn't exist
        stopifnot(
            # Stop if there are no images
            !is.null(SpatialExperiment::getImg(x)),
            # Stop if the image doesn't exist
            slice %in% SpatialExperiment::imgData(x)[1, "sample_id"],
            # Return error if there are no colnames in the object
            !is.null(colnames(x))
        )
        
        # If slice is null use the first slice
        img_df <- SpatialExperiment::imgData(x)
        if (is.null(slice))
            slice <- img_df[1, "sample_id"]
        
        # Scale factor to scale the coordinates
        sf <- img_df[img_df$sample_id == slice, "scaleFactor"]
        
        ## Extract spot barcodes
        barcodes <- colnames(x)
        
        ## Extract spatial coordinates
        # coord_df <- SpatialExperiment::spatialCoords(x)
        x <- as.matrix(SpatialExperiment::spatialCoords(x)[, c(1, 2)])
        
        ## Scale coordinates
        x <- x * sf
        
        ## Add barcodes to coord matrix & change colnames
        rownames(x) <- barcodes
    } else {
        stop("Couldn't extract image coordinates.
            Please check class(x) is SpatialExperiment, Seurat,
            dataframe or matrix")
    }
    return(x)
    
}
