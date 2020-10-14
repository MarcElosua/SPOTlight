#' This function takes in the H coefficient matrix object from and NMF object and returns plots to visualize the topic profiles between and within cell types
#'
#' @param h Object of class matrix; H coeficient matrix from NMF model.
#' @param train_cell_clust Object of class vector with cluster of the cells used to train the model.
#' @param clust_vr Object of class character. Name of the variable containing the cell clustering.
#' @return This function returns a list where the first element is a plot with the topic profiles of all the cell types and the 2nd element is a plot with the consensus topic profile per spot.
#' @export
#' @examples
#'

dot_plot_profiles_fun <- function(h,
                                  train_cell_clust,
                                  clust_vr) {

  suppressMessages(require(stringr)) # For the downsampling
  suppressMessages(require(dplyr))
  suppressMessages(require(tibble))
  suppressMessages(require(tidyr))
  suppressMessages(require(ggplot2))

  h_df <- data.frame(t(h))

  # Fix column names after converting to dataframe
  colnames(h_df) <- gsub(".", " ", colnames(h_df), fixed = TRUE)

  # Get proportions for each row
  h_ds <- round(h_df/rowSums(h_df), 4)
  h_ds[, clust_vr] <- train_cell_clust

  train_cells_plt <- h_ds %>%
    tibble::rowid_to_column("id") %>%
    tidyr::pivot_longer(cols = -c(all_of(clust_vr), id),
                        names_to = "topics",
                        values_to = "weights") %>%
    dplyr::group_by(!!! syms(clust_vr)) %>%
    dplyr::mutate(
      weights_txt = if_else(weights > 0.1, round(weights, 2), NULL)
    ) %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(aes(x = id, y = topics)) +
    ggplot2::geom_point(aes(size = weights, colour = weights)) +
    ggplot2::facet_wrap(as.formula(paste(clust_vr, "~ .")), scales = "free") +
    ggplot2::scale_color_continuous(low = "grey", high = "#59b371") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "NMF: Topic proportion within cell types") +
    ggplot2::theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.text = element_text(size = 15)) +
    ggplot2::scale_size(range = c(0, 5)) +
    ggplot2::guides(colour = guide_legend("Proportion"),
                    size = guide_legend("Proportion"))

  ct_topic_profiles <- h_ds %>%
    dplyr::group_by(!!! syms(clust_vr)) %>%
    dplyr::summarise_all(list(median)) %>%
    tibble::column_to_rownames(clust_vr)

  ct_topic_profiles <- ct_topic_profiles / rowSums(ct_topic_profiles)
  # In case a row is all 0
  ct_topic_profiles[is.na(ct_topic_profiles)] <- 0

  cell_type_plt <- round(ct_topic_profiles, 2) %>%
    tibble::rownames_to_column('Cell type') %>%
    tidyr::pivot_longer(cols = -`Cell type`, names_to = "Topics") %>%
    dplyr::mutate(
      value_txt = dplyr::if_else(value > 0.1, round(value, 2), NULL),
      Topics = factor(x = Topics,
                      levels = stringr::str_sort(colnames(ct_topic_profiles),
                                        numeric = TRUE))
    ) %>%
    ggplot2::ggplot(aes(x = `Cell type`, y = Topics)) +
    ggplot2::geom_point(aes(size = value, colour = value)) +
    ggplot2::scale_color_continuous(low = "grey", high = "#59b371") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "NMF: Topic profiles by cell type") +
    ggplot2::theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.text = element_text(size = 15)) +
    ggplot2::scale_size(range = c(0, 10)) +
    ggplot2::guides(colour = guide_legend("Proportion"),
                    size = guide_legend("Proportion"))

  return(list(train_cells_plt, cell_type_plt))
}

