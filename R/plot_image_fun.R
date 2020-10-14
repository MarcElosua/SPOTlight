#' This function takes in a seurat object, a path to an image and a slice name to load and plot.
#'
#' @param img_path: Object of class string containing the path pointing to the HE image.
#' @param img_alpha: Object of class numeric between 0-1 indicating the degree of transparency of the image
#' @return this function returns a plot of class gg
#' @export
#' @examples
#'

plot_image <- function(img_path,
                       img_alpha = 1) {

  # Check variables
  if (!is.character(img_path)) stop("ERROR: must be a character string!")
  if (!is.numeric(img_alpha)) stop("ERROR: img_alpha must be numeric between 0 and 1!")


  # Loading libraries
  suppressMessages(require(ggplot2))
  suppressMessages(require(imager))
  suppressMessages(require(dplyr))
  suppressMessages(require(tibble))

  ### Load histological image into R
  #### Extract file format, JPEG or PNG
  img_frmt <- tolower(stringr::str_sub(img_path, -4, -1))

  if(img_frmt %in% c(".jpg", "jpeg")) {
    img <- jpeg::readJPEG(img_path)
  } else if (img_frmt == ".png") {
    img <- png::readPNG(img_path)
  }

  # Convert image to grob object
  img_grob <- grid::rasterGrob(img,
                               interpolate = FALSE,
                               width = grid::unit(1, "npc"),
                               height = grid::unit(1, "npc"))

  ## Plot spatial scatterpie plot
  he_plt <- suppressMessages(ggplot2::ggplot() +
                                       ggplot2::annotation_custom(grob = img_grob,
                                                                  xmin = 0,
                                                                  xmax = 600,
                                                                  ymin = 0,
                                                                  ymax = -600) +
            ggplot2::scale_y_reverse() +
            ggplot2::ylim(600, 0) +
            ggplot2::xlim(0, 600) +
            cowplot::theme_half_open(11, rel_small = 1) +
            ggplot2::theme_void() +
            ggplot2::coord_fixed(ratio = 1,
                                 xlim = NULL,
                                 ylim = NULL,
                                 expand = TRUE,
                                 clip = "on"))

  return(he_plt)
}
