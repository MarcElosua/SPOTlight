#' This functions takes in a seurat object and returns the transposed count matrix with the highly variable genes selected
#'
#' @param se_obj Object of class Seurat with the data of interest
#' @return A dataframe with the ntop markers per cluster with the zscore of the logFC.
#' @export
#' @examples
#'

# se_obj: seurat object
# n: number of HVG desired, set by default to 3000
# noramlize: indicate if we want divie each cell by the sum of each row, in a way that the cells are comparable between them
#
# Returns:
# This function returns a sparse matrix object ready to pass to LDA function
######################################

prep_seobj_topic_fun <- function(se_obj){

  # Check variables
  if(is(se_obj)!="Seurat") {stop("ERROR: se_obj must be a Seurat object!")}

  #load required packages
  # cat("Loading packages...", sep="\n")
  suppressMessages(require(Seurat))
  suppressMessages(require(Matrix))

  # 1st get from the counts matrix from the RNA pocket, raw counts
  count_mtrx <- t(as.matrix(se_obj@assays$RNA@counts))

  # 2nd reconvert the matrix to sparse format again
  count_mtrx <- Matrix::Matrix(count_mtrx, sparse = T)

  return(count_mtrx)
}
