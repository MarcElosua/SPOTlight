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
  if (!"p_val" %in% colnames(markers)) stop("ERROR: p_val needs to be a variable in markers. markers must be the output of the FindAllMarkers Seurat function!")
  if (!"gene" %in% colnames(markers)) stop("ERROR: gene needs to be a variable in markers. markers must be the output of the FindAllMarkers Seurat function!")

  # load required packages
  suppressMessages(require(dplyr))

  # Check if we have logFC or avg_diff, if we have avg_diff assign it to avg_logFC to use it instead
  # if ("avg_diff" %in% colnames(markers) & ! "avg_logFC" %in% colnames(markers)) {
  #   markers$avg_logFC <- markers$avg_diff
  # }

  tmp_markers <- markers %>%
    dplyr::arrange(cluster, p_val) %>%
    dplyr::mutate(weight = 1 - p_val) %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(ntop) %>%
    dplyr::ungroup() %>%
    dplyr::select(gene, weight, p_val, cluster) %>%
    data.frame()

  return(tmp_markers)
}
