#' Run mixtures through the NMF model to get the cell type composition.
#'
#' @param nmf_mod Object of class dataframe obtained from the function Seurat::FindAllMarkers().
#' @param mixture_transcriptome Object of class matric of dimensions GENESxSPOTS
#' @param transf Object of class string indicatinf the transformation to normalize the count matrix: cpm (Counts per million), uv (unit variance), sct (Seurat::SCTransform), NULL (no transformation applied).
#' @param reference_profiles Object of class matrix containing the TOPICSxCELLS Coefficient matrix from where want to get the weights. It can be cell type profiles or cell specific profiles.
#' @return This function returns a matrix with the coefficients of the spatial mixtures.
#' @export
#' @examples
#'

mixture_deconvolution_nmf <- function(nmf_mod,
                                      mixture_transcriptome,
                                      transf,
                                      reference_profiles) {

  profile_mtrx <- predict_spatial_mixtures_nmf(nmf_mod = nmf_mod,
                               mixture_transcriptome = mixture_transcriptome,
                               transf = transf)

  # We add 1 extra column to add the residual error
  decon_mtrx <- matrix(data = NA, nrow = ncol(profile_mtrx), ncol = nrow(h) + 1)
  colnames(decon_mtrx) <- c(colnames(reference_profiles), "res_ss")

  for (i in seq_len(ncol(profile_mtrx))) {
    ## NNLS to get cell type composition
    nnls_pred <- nnls(A = reference_profiles, b = profile_mtrx[, i])

    ## get proportions of each cell type and multiply by 10 to get those cell types present > 10%, meaning representation of at least 1 cell
    comp <- nnls_pred$x / sum(nnls_pred$x) * 10
    comp[comp < 0.9] <- 0
    decon_mtrx[i, 1:(ncol(decon_mtrx) - 1)] <- comp
    decon_mtrx[i, ncol(decon_mtrx)] <- nnls_pred$deviance
  }

  res_ss_quant <- round(quantile(x = decon_mtrx[, ncol(decon_mtrx)],
                                 c(0.25, 0.5, 0.75),
                                 na.rm = TRUE), 5)

  print(sprintf("Quantiles of ss residuals are %s [%s - %s]",
                res_ss_quant[2], res_ss_quant[1], res_ss_quant[3]))

  return(decon_mtrx)
}
