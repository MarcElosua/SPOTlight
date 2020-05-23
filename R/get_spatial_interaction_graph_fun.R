#' This function takes in deconvolution matrix and returns a network igraph object
#'
#' @param decon_mtrx Object of class matrix, SPOTxCELL-TYPE, with the cell-type proportions per spot.
#' @export
#' @examples
#'


get_spatial_interaction_graph <- function(decon_mtrx) {

  # Check variables
  if (!is.matrix(decon_mtrx)) stop("ERROR: decon_mtrx must be a matrix object!")

  # Require needed libraries
  suppressMessages(require(igraph))

  if (is.null(colnames(decon_mtrx))) {
    colnames(decon_mtrx) <- as.character(1:ncol(decon_mtrx))
  }
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
  ntwrk_mtrx <- ntwrk_mtrx[ntwrk_mtrx[, 3] != "0", ]
  # add column with scaled values
  ntwrk_mtrx <- cbind(ntwrk_mtrx, scale(as.numeric(ntwrk_mtrx[, 3]), center = 1))


  # data <- matrix(sample(0:1, 100, replace=TRUE, prob=c(0.8,0.2)), nc=10)
  links <- data.frame(
    source = ntwrk_mtrx[, 1],
    target = ntwrk_mtrx[, 2],
    importance = as.numeric(ntwrk_mtrx[, 4])
  )
  nodes <- data.frame(name=colnames(decon_mtrx))

  network <- igraph::graph_from_data_frame(d = links,
                                           vertices = nodes,
                                           directed = FALSE)

  return(network)
}
