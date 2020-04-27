#' This function takes in the H coefficient matrix object from and NMF object and returns plots to visualize the topic profiles between and within cell types
#'
#' @param h Object of class matrix; H coeficient matrix from NMF model.
#' @param train_cell_clust Object of class vector with cluster of the cells used to train the model.
#' @param clust_vr Object of class character. Name of the variable containing the cell clustering.
#' @return This function returns a list where the first element is a matrix with the topic profiles of all possible combinations and the 2nd element is the cell composition of each spot.
#' @export
#' @examples
#'

dot_plot_profiles_fun <- function(h,
                                  train_cell_clust,
                                  clust_vr) {

  h_ds <- data.frame(t(h))
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
    ggplot(aes(x = id, y = topics)) +
    geom_point(aes(size = weights, colour = weights)) +
    facet_wrap(as.formula(paste(clust_vr, "~ .")), scales = "free") +
    scale_color_continuous(low = "grey", high = "Green") +
    theme_classic() +
    labs(title = "NMF: Topic proportion within cell types") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text.x = element_text(angle = 45, vjust = 0.5),
      axis.text = element_text(size = 15)) +
    scale_size(range = c(0, 20)) +
    guides(colour = guide_legend(""), size = guide_legend(""))

  ct_topic_profiles <- h_ds %>%
    dplyr::group_by(!!! syms(clust_vr)) %>%
    dplyr::summarise_all(list(median)) %>%
    tibble::column_to_rownames(clust_vr) %>%
    as.matrix()

  cell_type_plt <- round(ct_topic_profiles, 2) %>%
    data.frame() %>%
    tibble::rownames_to_column('Cell type') %>%
    tidyr::pivot_longer(cols = -`Cell type`, names_to = "Topics") %>%
    mutate(
      value_txt = if_else(value > 0.1, round(value, 2), NULL),
      Topics = factor(x = Topics,
                      levels = str_sort(colnames(ct_topic_profiles),
                                        numeric = TRUE))
    ) %>%
    ggplot(aes(x = `Cell type`, y = Topics)) +
    geom_point(aes(size = value, colour = value)) +
    scale_color_continuous(low = "grey", high = "Green") +
    theme_classic() +
    labs(title = "NMF: Topic profiles by cell type") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text.x = element_text(angle = 45, vjust = 0.5),
      axis.text = element_text(size = 15)) +
    scale_size(range = c(0, 15)) +
    guides(colour = guide_legend(""), size = guide_legend(""))

  return(list(train_cells_plt, cell_type_plt))
}

