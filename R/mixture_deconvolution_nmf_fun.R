#' Run mixtures through the NMF model to get the cell type composition.
#'
#' @param nmf_mod Object of class NMF containing the trained nNMF model.
#' @param mixture_transcriptome Object of class matric of dimensions GENESxSPOTS
#' @param transf Transformation to normalize the count matrix: uv (unit variance), raw (no transformation applied). By default UV.
#' @param reference_profiles Object of class matrix containing the TOPICSxCELLS Coefficient matrix from where want to get the weights. It can be cell type profiles or cell specific profiles.
#' @param min_cont Object of class numeric; Indicates the minimum contribution we expect from a cell in that spot. Since we're working with proportions by setting 0.01, by default, means that we will accept those cell types whose weight coefficient is at least 1\% of the total.
#' @return This function returns a matrix with the coefficients of the spatial mixtures.
#' @export
#' @examples
#'

mixture_deconvolution_nmf <- function(nmf_mod,
                                      mixture_transcriptome,
                                      transf,
                                      reference_profiles,
                                      min_cont = 0.01) {

  # Check variables
  if (!is(nmf_mod, "NMF")) stop("ERROR: nmf_mod must be an NMF object!")
  if (!is.character(transf)) stop("ERROR: transf must be a character string!")
  if (!is.matrix(reference_profiles)) stop("ERROR: reference_profiles must be a matrix!")
  if (!is.numeric(min_cont)) stop("ERROR: min_cont must be numeric!")

  # Loading libraries
  suppressMessages(require(nnls))

  profile_mtrx <- predict_spatial_mixtures_nmf(nmf_mod = nmf_mod,
                               mixture_transcriptome = mixture_transcriptome,
                               transf = transf)

  # We add 1 extra column to add the residual error
  decon_mtrx <- matrix(data = NA,
                       nrow = ncol(profile_mtrx),
                       ncol = ncol(reference_profiles) + 1)
  colnames(decon_mtrx) <- c(colnames(reference_profiles), "res_ss")

  # create progress bar
  print("Deconvoluting spots")
  total <- ncol(profile_mtrx)
  pb <- txtProgressBar(min = 0, max = total, style = 3)

  for (i in seq_len(ncol(profile_mtrx))) {
    ## NNLS to get cell type composition
    nnls_pred <- nnls::nnls(A = reference_profiles, b = profile_mtrx[, i])
    weights <- nnls_pred$x

    ## get proportions of each cell type
    comp <- weights / sum(weights)

    ## Remove cell types not contributing the minimum
    comp[comp < min_cont] <- 0
    weights[comp < min_cont] <- 0

    ### Updated proportions after filtering out minimum contributions
    comp_prop <- comp / sum(comp)
    comp_prop[is.na(comp_prop)] <- 0

    ## Get Total sum of squares
    fit_null <- 0
    tot_ss <- sum((profile_mtrx[, i] - fit_null) ^ 2)

    ## Get % of unexplained residuals
    unexpl_ss <- nnls_pred$deviance / tot_ss

    decon_mtrx[i, 1:(ncol(decon_mtrx) - 1)] <- comp_prop
    decon_mtrx[i, ncol(decon_mtrx)] <- unexpl_ss

    # update progress bar
    setTxtProgressBar(pb, i)
  }
  # Close progress bar
  close(pb)

  return(decon_mtrx)
}
