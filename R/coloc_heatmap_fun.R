#' This function takes in deconvolution matrix and returns a network igraph object
#'
#' @param decon_mtrx Object of class matrix, SPOTxCELL-TYPE, with the cell-type proportions per spot.
#' @return This function returns a heatmap representing how the cell type colocalize on the tissue
#' @export
#' @examples
#'


get_colocalization_heatmap <- function(decon_mtrx) {
  comb_id <- arrangements::combinations(x = colnames(decon_mtrx), k = 2, replace = F)
  comb_id_str <- paste(comb_id[, 1], comb_id[, 2], sep = "_")
  comb_val <- matrix(data = 0, nrow = nrow(comb_id), ncol = 1)
  rownames(comb_val) <- comb_id_str

  # Iterate over all the spot's predicted cell composition
  for (i in seq_len(nrow(decon_mtrx))) {
    mtp_row <- decon_mtrx[i, ]
    mtp_row_sub <- mtp_row[mtp_row != 0]

    # If there is only one cell type ignore the loop
    if (length(names(mtp_row_sub)) > 1) {

      # Iterate over all the cell types within that spot and fill the comb_val matrix
      # names(mtp_row_sub)[(pos+1):length(mtp_row_sub)] - Set iterator pos to avoid counting things twice
      # ii in names(mtp_row_sub)[-length(mtp_row_sub)] - Don't iterate over the last one since it will have been already counted by all the previous ones
      pos <- 1
      for (ii in names(mtp_row_sub)[-length(mtp_row_sub)]) {

        for(iii in names(mtp_row_sub)[(pos+1):length(mtp_row_sub)]){
          tmp_id <- paste(ii, iii, sep = "_")
          comb_val[rownames(comb_val) == tmp_id,] = comb_val[rownames(comb_val) == tmp_id,] + 1
        }
        pos = pos + 1
      }
    }
  }
  # Join matrices and scale comb_val centering it around 1
  ntwrk_mtrx <- cbind(comb_id, comb_val)
  # Remove rows belonging to cell types not interacting
  # ntwrk_mtrx <- ntwrk_mtrx[ntwrk_mtrx[, 3] != "0", ]
  # add column with scaled values
  # ntwrk_mtrx <- cbind(ntwrk_mtrx, scale(as.numeric(ntwrk_mtrx[, 3]), center = 1))

  tmp_df <- data.frame(ntwrk_mtrx)
  colnames(tmp_df) <- c("ct1", "ct2", "colocalizations")
  # Get on how many spots we have each cell type
  total_spots <- data.frame(total_spots = colSums(decon_mtrx > 0)) %>%
    tibble::rownames_to_column("cell_type")


  tmp_df <- tmp_df %>%
    dplyr::left_join(total_spots, by = c("ct1" = "cell_type")) %>%
    dplyr::left_join(total_spots, by = c("ct2" = "cell_type"), suffix = c("_1", "_2")) %>%
    mutate(
      colocalizations = as.numeric(colocalizations),
      coloc_prop_1 = colocalizations / total_spots_1,
      coloc_prop_2 = colocalizations / total_spots_2,
      ) %>%
    tidyr::pivot_longer(cols = c("coloc_prop_1", "coloc_prop_2"),
                        names_to = "ct_denom",
                        values_to = "coloc_prop") %>%
    mutate(coloc_prop = if_else(is.na(coloc_prop), 0, coloc_prop))

  tmp_df1 <- tmp_df %>% dplyr::filter(ct_denom == "coloc_prop_1")
  coloc_df <-  tmp_df %>% dplyr::filter(ct_denom == "coloc_prop_2") %>%
    dplyr::mutate(ct3 = ct1,
           ct1 = ct2,
           ct2 = ct3) %>%
    dplyr::select(-ct3) %>%
    rbind(tmp_df1)


  hm_plt <- ggplot(coloc_df,
         aes(x = ct1,
             y = ct2,
             fill = coloc_prop,
             alpha = coloc_prop)) +
  geom_tile() +
  labs(title = "Colocalization heatmap",
       fill = "Contacts",
       alpha = "Contacts",
       x = "Reference cell type",
       y = "Cell type colocalizing with") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
  ) +
  scale_fill_gradient(low = "lightgrey", high = "red", guide = "legend") +
  scale_alpha_continuous(guide = "legend")

  return(hm_plt)
}
