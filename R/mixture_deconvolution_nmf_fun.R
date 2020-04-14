#' Run mixtures through the NMF model to get the cell type composition.
#'
#' @param nmf_mod Object of class dataframe obtained from the function Seurat::FindAllMarkers().
#' @param mixture_transcriptome Object of class matric of dimensions GENESxSPOTS
#' @param transf Object of class string indicatinf the transformation to normalize the count matrix: cpm (Counts per million), uv (unit variance), raw (no transformation applied).
#' @param reference_profiles Object of class matrix containing the TOPICSxCELLS Coefficient matrix from where want to get the weights. It can be cell type profiles or cell specific profiles.
#' @param min_cont Object of class numeric; Indicates the minimum contribution we expect from a cell in that spot. Since we're working with proportions by setting 0.09, by default, means that we will accept those cell types whose weight coefficient is at least 9% of the total.
#' @return This function returns a matrix with the coefficients of the spatial mixtures.
#' @export
#' @examples
#'

mixture_deconvolution_nmf <- function(nmf_mod,
                                      mixture_transcriptome,
                                      transf,
                                      reference_profiles,
                                      min_cont = 0.09) {

  # Check variables
  if (!is(h, "matrix")) stop("ERROR: h must be a matric object!")
  if (! is(train_cell_clust, "vector")) stop("ERROR: train_cell_clust must be a vector/list object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")

  # Loading libraries
  suppressMessages(require(nnls))

  profile_mtrx <- predict_spatial_mixtures_nmf(nmf_mod = nmf_mod,
                               mixture_transcriptome = mixture_transcriptome,
                               transf = transf)

  # We add 1 extra column to add the residual error
  decon_mtrx <- matrix(data = NA, nrow = ncol(profile_mtrx), ncol = ncol(reference_profiles) + 1)
  colnames(decon_mtrx) <- c(colnames(reference_profiles), "res_ss")

  # create progress bar
  print("Deconvoluting spots")
  total <- ncol(profile_mtrx)
  pb <- txtProgressBar(min = 0, max = total, style = 3)

  for (i in seq_len(ncol(profile_mtrx))) {
    ## NNLS to get cell type composition
    nnls_pred <- nnls::nnls(A = reference_profiles, b = profile_mtrx[, i])

    ## get proportions of each cell type and multiply by 10 to get those cell types present > 10%, meaning representation of at least 1 cell
    comp <- nnls_pred$x / sum(nnls_pred$x)
    comp[comp < min_cont] <- 0
    decon_mtrx[i, 1:(ncol(decon_mtrx) - 1)] <- comp*10
    decon_mtrx[i, ncol(decon_mtrx)] <- nnls_pred$deviance

    # update progress bar
    setTxtProgressBar(pb, i)
  }
  # Close progress bar
  close(pb)

  return(decon_mtrx)
}
