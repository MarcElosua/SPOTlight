#' This function takes in a seurat object and cell types of interest and returns a scatterpie plot with each spot situated in its spatial location.
#'
#' @param se_obj: Object of class Seurat with the spatial data and the cell type proportions in the metadata.
#' @param cell_types_all: Object of class vector containing the names of all the cell types.
#' @param img_path: Object of class string pointing to the HE image.
#' @param cell_types_interest: Object of class vector containing the cell types of interest you want to plot. By setting this parameters only spots containing at least one of these cell types will be plotted. By default, NULL, it will be assumed to be the same as cell_types_all.
#' @param return_legend: Object of class logical. By default (FALSE) it will return the legend along with the plot. If TRUE it will return the plot and the legend separately as grob object, this is a useful option if you want to do image compositions later on.
#' @param slice: Object of class character, name of the slice image to load as found in se_obj@images, by default it will grab the first one on the list.
#' @param img_alpha: Object of class numeric between 0-1 indicating the degree of transparency of the image.
#' @param scatterpie_alpha: Object of class numeric between 0-1 indicating the degree of transparency of the scatterpie.
#' @param pie_scale: Object of class numeric containing the size of the pie charts.
#' @param col_df: object of class dataframe containing a color for each cell type, the first column is the cell type and the second is the color assigned to it.
#' @return this function returns a plot composition of class gg
#' @export
#' @examples

spatial_scatterpie <- function(se_obj,
                               cell_types_all,
                               img_path,
                               cell_types_interest = NULL,
                               return_legend = FALSE,
                               slice = NULL,
                               img_alpha = 1,
                               scatterpie_alpha = 1,
                               pie_scale = 0.4,
                               col_df = NULL) {

  # Check variables
  if (!is(se_obj, "Seurat")) stop("ERROR: se_obj must be a Seurat object!")
  if (! is(cell_types_all, "vector")) stop("ERROR: cell_types_all must be a vector/list object!")
  if (!is.character(img_path)) stop("ERROR: must be a character string!")
  if (!(is(cell_types_interest, "vector") | is.null(cell_types_interest))) stop("ERROR: cell_types_interest must be a vector/list object or NULL!")
  if (!is(return_legend, "logical")) stop("ERROR: return_legend must be logical!")
  if (is.character(slice) & !slice %in% names(se_obj@images)) stop("ERROR: slice is not in names(se_obj@images)!")
  if (!is.numeric(img_alpha)) stop("ERROR: img_alpha must be numeric between 0 and 1!")
  if (!is.numeric(scatterpie_alpha)) stop("ERROR: scatterpie_alpha must be numeric between 0 and 1!")
  if (!is.numeric(pie_scale)) stop("ERROR: pie_scale must be numeric between 0 and 1!")
  if (!(is.data.frame(col_df) | is.null(col_df))) stop("ERROR: col_df must be a dataframe or NULL!")

  # Loading libraries
  suppressMessages(require(ggplot2))
  suppressMessages(require(cowplot))
  suppressMessages(require(imager))
  suppressMessages(require(dplyr))
  suppressMessages(require(tibble))

  # Load the image, plot the scatterplot and overlap them:
  # Load the image + plot the image
  he_plt <- plot_image(img_path = img_path,
                       img_alpha = img_alpha)

  # Plot the scatterplot
  scatterpie_plt <- scatterpie_plot(se_obj = se_obj,
                                    cell_types_all = cell_types_all,
                                    slice = slice,
                                    scatterpie_alpha = scatterpie_alpha,
                                    cell_types_interest = cell_types_interest,
                                    col_df = col_df) +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

  if (return_legend) {
    legend_grob <- get_legend(scatterpie_plt)
    # Align and superopose plots
    scatterpie_plt <- scatterpie_plt + theme(legend.position="none")
  }
  # Align and superopose plots
  aligned_plts <- align_plots(he_plt,
                              scatterpie_plt,
                              align = "hv",
                              axis = "tblr")

  spatial_superposed <- ggdraw(aligned_plts[[1]]) + draw_plot(aligned_plts[[2]])

  if (return_legend) {
    return(list(spatial_superposed, legend_grob))
  } else {
    return(spatial_superposed)
  }

}
