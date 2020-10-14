#' This functions takes in a seurat object and returns the transposed count matrix with the highly variable genes selected
#'
#' @param se_obj Object of class Seurat with the data of interest
#' @return This function returns a sparse matrix object.
#' @export
#' @examples
#'

prep_seobj_topic_fun <- function(se_obj) {

  # Check variables
  if (is(se_obj) != "Seurat") stop("ERROR: se_obj must be a Seurat object!")

  #load required packages
  suppressMessages(require(Seurat))
  suppressMessages(require(Matrix))

  # 1st get from the counts matrix from the RNA pocket, raw counts
  count_mtrx <- t(as.matrix(se_obj@assays$RNA@counts))

  # 2nd reconvert the matrix to sparse format again
  count_mtrx <- Matrix::Matrix(count_mtrx, sparse = T)

  return(count_mtrx)
}
