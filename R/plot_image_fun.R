#' This function takes in a seurat object, a path to an image and a slice name to load and plot.
#'
#' @param img_path: Object of class string pointing to the HE image.
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

  spatial_img <- imager::load.image(file = img_path)

  ### Plot image ###
  he_plt <- spatial_img %>%
    as.data.frame(wide = "c") %>%
    mutate(rgb.val = rgb(c.1, c.2, c.3)) %>%
    ggplot(., aes(x, y)) +
      geom_raster(aes(fill = rgb.val), alpha = img_alpha) +
      scale_fill_identity() +
      scale_y_reverse() +
      theme_void() +
      theme(plot.margin = unit(c(0, 0, 0, 0), "line")) +
      coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

  return(he_plt)
}
