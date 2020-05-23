#' This function takes in a seurat object and cell types of interest and returns a scatterpie plot with each spot situated in its spatial location.
#'
#' @param se_obj: Object of class Seurat with the spatial data and the cell type proportions in the metadata.
#' @param cell_types_all: Object of class vector containing the names of all the cell types.
#' @param slice: Object of class character, name of the slice image to load as found in se_obj@images, by default it will grab the first one on the list.'
#' @param scatterpie_alpha: Object of class numeric between 0-1 indicating the degree of transparency of the scatterpie
#' @param cell_types_interest: Object of class vector containing the cell types of interest you want to plot. By setting this parameters only spots containing at least one of these cell types will be plotted. By default, NULL, it will be assumed to be the same as cell_types_all.
#' @param pie_scale: Object of class numeric containing the size of the pie charts.
#' @param col_df: object of class dataframe containing a color for each cell type, the first column is the cell type and the second is the color assigned to it.
#' @return this function returns a plot composition of class gg
#' @export
#' @examples
#'

scatterpie_plot <- function(se_obj,
                            cell_types_all,
                            slice = NULL,
                            scatterpie_alpha = 1,
                            cell_types_interest = NULL,
                            pie_scale = 0.4,
                            col_df = NULL) {

  # Check variables
  if (!is(se_obj, "Seurat")) stop("ERROR: se_obj must be a Seurat object!")
  if (! is(cell_types_all, "vector")) stop("ERROR: cell_types_all must be a vector/list object!")
  if (!(is(cell_types_interest, "vector") | is.null(cell_types_interest))) stop("ERROR: cell_types_interest must be a vector/list object or NULL!")
  # if (is.null(slice) | (!is.null(slice) && !slice %in% names(se_obj@images))) warning("Warning: slice is not in names(se_obj@images)!, 1st name will be used!")
  if (!is.numeric(scatterpie_alpha)) stop("ERROR: scatterpie_alpha must be numeric between 0 and 1!")
  if (!is.numeric(pie_scale)) stop("ERROR: pie_scale must be numeric between 0 and 1!")
  if (!(is.data.frame(col_df) | is.null(col_df))) stop("ERROR: col_df must be a dataframe or NULL!")

  # Loading libraries
  suppressMessages(require(ggplot2))
  suppressMessages(require(cowplot))
  suppressMessages(require(imager))
  suppressMessages(require(dplyr))
  suppressMessages(require(tibble))

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

  # Plot the scatterplot
  scatterpie_plt <- suppressMessages(ggplot() +
                     scatterpie::geom_scatterpie(data = spatial_coord,
                                                 aes(x = imagecol_scaled,
                                                     y = imagerow_scaled),
                                                 cols = cell_types_all,
                                                 color = NA,
                                                 alpha = scatterpie_alpha,
                                                 pie_scale = pie_scale) +
                     scale_y_reverse() +
                     ylim(600, 0) +
                     xlim(0, 600) +
                     theme_half_open(11, rel_small = 1) +
                     theme_void())

  if (!is.null(col_df)) {
    # Select which cell types are going to be plotted
    ct_keep <- sort(names(which(colSums(spatial_coord[, cell_types_all]) > 0)))

    # Subset colours to use
    col_df[, 1] <- gsub(pattern = "[[:punct:]]|[[:blank:]]",
                        replacement = ".",
                        x = col_df[, 1],
                        perl = TRUE)

    col_vec <- col_df[col_df[, 1] %in% ct_keep, 2]

    # Add colours to plot
    scatterpie_plt <- scatterpie_plt + scale_fill_manual(values = col_vec)
  }

  return(scatterpie_plt)

}
