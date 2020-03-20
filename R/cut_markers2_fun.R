#' Select the top n cluster markers by ranking the genes by pvalue and then by avglogFC
#'
#' @param markers A data.frame as in the output of the FindAllMarkers Seurat function.
#' @param ntop The number of top markers you want.
#' @return A dataframe with the ntop markers per cluster with the zscore of the logFC.
#' @export
#' @examples
#'

cut_markers2 <- function(markers,
                         ntop) {

  # Check variables
  if (!is.data.frame(markers)) stop("ERROR: markers must be a data.frame object!")
  if (!is.numeric(ntop)) stop("ERROR: ntop must be a an integer!")
  if (!"cluster" %in% colnames(markers)) stop("ERROR: cluster needs to be a variable in markers. markers must be the output of the FindAllMarkers Seurat function!")
  if (!"p_val_adj" %in% colnames(markers)) stop("ERROR: p_val_adj needs to be a variable in markers. markers must be the output of the FindAllMarkers Seurat function!")
  if ((! "avg_logFC" %in% colnames(markers)) & (! "avg_diff" %in% colnames(markers))) stop("ERROR: avg_logFC or avg_diff needs to be a variable in markers. markers must be the output of the FindAllMarkers Seurat function!")
  if (!"gene" %in% colnames(markers)) stop("ERROR: gene needs to be a variable in markers. markers must be the output of the FindAllMarkers Seurat function!")

  # load required packages
  suppressMessages(require(dplyr))

  # Check if we have logFC or avg_diff, if we have avg_diff assign it to avg_logFC to use it instead
  if ("avg_diff" %in% colnames(markers) & ! "avg_logFC" %in% colnames(markers)) {
    markers$avg_logFC <- markers$avg_diff
  }

  tmp_markers <- markers %>%
    dplyr::arrange(cluster, p_val_adj, avg_logFC) %>%
    dplyr::mutate(
      # If there are + or - infinite values set them as the highest value or lowest value, we are just ranking them so the absolute value isn't crucial
      avg_logFC_mod = dplyr::if_else(is.finite(avg_logFC), avg_logFC, NA_real_),
      avg_logFC = dplyr::if_else(is.finite(avg_logFC), avg_logFC,
                                 dplyr::if_else(avg_logFC > 0,
                                                max(avg_logFC_mod, na.rm = TRUE) + 5,
                                                min(avg_logFC_mod, na.rm = TRUE) - 5)),
      # Scale the logFC, we will use this number to seed the model
      # We are scaling by dividing all the avg_logFC by the max
      logFC_z = scale(x = avg_logFC, center = F, scale = max(avg_logFC))
      ) %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(ntop) %>%
    dplyr::ungroup() %>%
    dplyr::select(gene, logFC_z, cluster) %>%
    data.frame()

  return(tmp_markers)
}
