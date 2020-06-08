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
#'

spatial_scatterpie_alt <- function(se_obj,
                               cell_types_all,
                               img_path,
                               cell_types_interest = NULL,
                               return_legend = FALSE,
                               slice = NULL,
                               img_alpha = 1,
                               scatterpie_alpha = 1,
                               pie_scale = 1,
                               col_df = NULL) {

  # Check variables
  if (!is(se_obj, "Seurat")) stop("ERROR: se_obj must be a Seurat object!")
  if (! is(cell_types_all, "vector")) stop("ERROR: cell_types_all must be a vector/list object!")
  if (!is.character(img_path)) stop("ERROR: must be a character string!")
  if (!(is(cell_types_interest, "vector") | is.null(cell_types_interest))) stop("ERROR: cell_types_interest must be a vector/list object or NULL!")
  if (!is(return_legend, "logical")) stop("ERROR: return_legend must be logical!")
  # if (is.character(slice) & !slice %in% names(se_obj@images)) stop("ERROR: slice is not in names(se_obj@images)!")
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
  suppressMessages(require(Spaniel))

  metadata_ds <- data.frame(se_obj@meta.data)

  ## Change column names for consistency ##
  # [[:punct:]] - Any punctuation character: ! ' # S % & ' ( ) * + , - . / : ; < = > ? @ [ / ] ^ _ { | } ~
  colnames(metadata_ds) <- gsub(pattern = "[[:punct:]]|[[:blank:]]",
                                replacement = ".",
                                x = colnames(metadata_ds),
                                perl = TRUE)

  cell_types_all <- gsub(pattern = "[[:punct:]]|[[:blank:]]",
                         replacement = ".",
                         x = cell_types_all,
                         perl = TRUE)

  if (is.null(cell_types_interest)) {
    cell_types_interest <- cell_types_all
  } else {
    cell_types_interest <- gsub(pattern = "[[:punct:]]|[[:blank:]]",
                                replacement = ".",
                                x = cell_types_interest,
                                perl = TRUE)
  }

  # If not all cell types are in the cell types of interest we only want to keep those spots which have at least one of the cell types of interest
  if (!all(cell_types_all %in% cell_types_interest)) {

    metadata_ds <- metadata_ds %>%
      tibble::rownames_to_column("ID") %>%
      dplyr::mutate(rsum = rowSums(.[, cell_types_interest, drop = FALSE])) %>%
      dplyr::filter(rsum != 0) %>%
      dplyr::select("ID") %>%
      dplyr::left_join(metadata_ds %>% rownames_to_column("ID"),by = "ID") %>%
      tibble::column_to_rownames("ID")
  }

  ## If slice is not selected set it to the first element in the list of slices
  if (is.null(slice) | (!is.null(slice) && !slice %in% names(se_obj@images))) {
    slice <- names(se_obj@images)[1]
    print(sprintf("Using slice %s", slice))
  }

  ## Preprocess data
  spatial_coord <- data.frame(se_obj@images[[slice]]@coordinates) %>%
    tibble::rownames_to_column("ID") %>%
    dplyr::mutate(imagerow_scaled = imagerow * se_obj@images[[slice]]@scale.factors$lowres,
                  imagecol_scaled = imagecol * se_obj@images[[slice]]@scale.factors$lowres) %>%
    dplyr::inner_join(metadata_ds %>% tibble::rownames_to_column("ID"), by = "ID")

  ### Load histological image into R
  img <- png::readPNG(imgFile)
  img_grob <- grid::rasterGrob(img,
                        interpolate = FALSE,
                        width = grid::unit(1, "npc"),
                        height = grid::unit(1, "npc"))

  ## Plot spatial scatterpie plot
  scatterpie_plt <- suppressMessages(ggplot2::ggplot() +
    ggplot2::annotation_custom(grob = img_grob,
                               xmin = 0,
                               xmax = 600,
                               ymin = 0,
                               ymax = -600) +
    scatterpie::geom_scatterpie(data = spatial_coord,
                                aes(x = imagecol_scaled,
                                    y = imagerow_scaled),
                                cols = cell_types_all,
                                color = NA,
                                alpha = scatterpie_alpha,
                                pie_scale = pie_scale) +
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



}
