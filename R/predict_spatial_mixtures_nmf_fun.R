#' Run test spots through the basis to get the pertinent coefficients. To do this for every spot we are going to set up a system of linear equations where we need to find the coefficient, we will use non-negative least squares to determine the best coefficient fit. Returns a matrix of KxnSPOTS dimensions.
#'
#' @param nmf_mod Object of class dataframe obtained from the function Seurat::FindAllMarkers().
#' @param mixture_transcriptome Object of class matric of dimensions GENESxSPOTS
#' @param transf Transformation to normalize the count matrix: cpm (Counts per million), uv (unit variance), sct (Seurat::SCTransform), NULL (no transformation applied).
#' @return This function returns a matrix with the coefficients of the spatial mixtures.
#' @export
#' @examples
#'

predict_spatial_mixtures_nmf <- function(nmf_mod,
                                         mixture_transcriptome,
                                         transf) {

  ## Extract genes used in w, if there are genes not present add them with all 0
  keep_genes <- rownames(basis(nmf_mod))[rownames(basis(nmf_mod)) %in% rownames(mixture_transcriptome)]
  mixture_transcriptome_subs <- as.matrix(mixture_transcriptome[keep_genes, ])

  if (transf == "cpm") {
    count_mtrx <- edgeR::cpm(mixture_transcriptome,
                             normalized.lib.sizes = FALSE)

  } else if (transf == "uv") {
    count_mtrx <- scale(t(mixture_transcriptome),
                        center = FALSE,
                        scale = apply(t(counts), 2, sd, na.rm = TRUE))
    count_mtrx <- t(count_mtrx)

  } else if (transf == "sct") {
    # Can't use scale.data since it has negative values
    count_mtrx <- mixture_transcriptome

  } else if (is.null(transf)) {
    count_mtrx <- mixture_transcriptome

  }

  ##### Extract Basis matrix W #####
  W <- basis(nmf_mod)
  H <- coef(nmf_mod)

  coef_pred <- matrix(data = NA,
                    nrow = ncol(W),
                    ncol = ncol(count_mtrx))
  colnames(coef_pred) <- colnames(count_mtrx)

  ##### Perform NNLS to get coefficients #####
  for (i in seq_len(ncol(count_mtrx))) {
    nnls_pred <- nnls(A = W, b = count_mtrx[, i])
    coef_pred[, i] <- nnls_pred$x
  }



  return(coef_pred)
}
