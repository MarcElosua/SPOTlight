#' Select the top n cluster markers by ranking the genes by pvalue and then by avglogFC
#'
#' @param markers A data.frame as in the output of the FindAllMarkers Seurat function.
#' @param ntop The number of top markers you want.
#' @return A dataframe with the ntop markers per cluster with the zscore of the logFC.
#' @export
#' @examples
#' 

cut_markers2 <- function(markers, ntop){
  
  # Check variables
  if(!is.data.frame(markers)) {stop("ERROR: markers must be a data.frame object!")}
  if(!is.integer(ntop)) {stop("ERROR: ntop must be a an integer!")}
  if("cluster" %in% colnames(markers)) {stop("ERROR: cluster needs to be a variable in markers. markers must be the output of the FindAllMarkers Seurat function!")}
  if("p_val_adj" %in% colnames(markers)) {stop("ERROR: p_val_adj needs to be a variable in markers. markers must be the output of the FindAllMarkers Seurat function!")}
  if("avg_logFC" %in% colnames(markers)) {stop("ERROR: avg_logFC needs to be a variable in markers. markers must be the output of the FindAllMarkers Seurat function!")}
  if("gene" %in% colnames(markers)) {stop("ERROR: gene needs to be a variable in markers. markers must be the output of the FindAllMarkers Seurat function!")}
  
  #load required packages
  cat("Loading packages...", sep="\n")
  suppressMessages(require(dplyr))

  tmp_markers <- markers %>%
    filter(p_val_adj < 0.01) %>%
    arrange(cluster,p_val_adj,avg_logFC) %>%
    mutate(logFC_z = abs(scale(avg_logFC))) %>% 
    group_by(cluster) %>% 
    top_n(ntop) %>%
    ungroup() %>% 
    dplyr::select(gene, logFC_z, cluster) %>% 
    data.frame()
  
  return(tmp_markers)
}
