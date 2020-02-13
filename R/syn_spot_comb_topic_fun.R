#' This function takes in the cluster profiles and returns a sotchastic combination of all the possible ones
#'
#' @param lda_mod Object of class LDA_Gibbs.
#' @param se_obj Object of class Seurat.
#' @param clust_vr Object of class character. Name of the variable containing the cell clustering.
#' @param verbose Object of class Logical determining if progress should be reported or not (TRUE by default).
#' @return This function returns a list where the first element is a matric with the topic profiles of all possible combinations and the 2nd element is the cell composition of each spot.
#' @export
#' @examples
#'

syn_spot_comb_topic <- function(lda_mod, se_obj, clust_vr, verbose = TRUE) {

  # Check variables
  if (is(lda_mod)[[1]] != "LDA_Gibbs") stop("ERROR: lda_mod must be a LDA_Gibbs object!")
  if (is(se_obj) != "Seurat") stop("ERROR: se_obj must be a Seurat object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  if (!is.logical(verbose)) stop("ERROR: verbose must be a logical object!")

  # Load needed packages
  suppressMessages(require(Matrix))
  suppressMessages(require(arrangements))
  suppressMessages(require(progress))
  suppressMessages(require(purrr))
  suppressMessages(require(topicmodels))


  #### Calculate topic profiles for every cluster ####
  clust_profiles <- topic_profile_per_cluster(lda_mod = lda_mod,
                                              se_obj = se_obj,
                                              clust_vr = clust_vr)
  round(clust_profiles, 4)

  if (class(clust_profiles) != "matrix") clust_profiles <- as.matrix(clust_profiles)

  # If a cluster is 0 change it to 1
  if (sum(grepl(pattern = "0", rownames(clust_profiles))) != 0) {
    rownames(clust_profiles) <- as.character(as.numeric(rownames(clust_profiles)) + 1)
  }

  # Compute all possible combinations
  k_sub <- 8

  # First compute how many combinations are possible
  total_comb <- arrangements::ncombinations(x = c(0:nrow(clust_profiles)), k = k_sub, replace = TRUE)

  # If the total number of combinations exceeds 500K with k_sub = 8 we will select 1M random combinations
  # 18 cluster with k_sub 8 gives ~1M possible combinations,
  if (total_comb > 1e6) {
    # select 1M from all the k_sub = 8 possible comb and the rest will be discarded.
    if (verbose) print(sprintf("The total number of possible combinations selecting up to 8 cells per spot is: %s, using %s random combinations", total_comb, 1e6))
    comb <- arrangements::combinations(x = c(0 : nrow(clust_profiles)),
                                       k = k_sub, replace = TRUE, nsample = 1e6)
  } else {
    # Do all possible combinatinos
    if (verbose) print(sprintf("Generating all synthetic spot combinations: %s", total_comb))
    comb <- arrangements::combinations(x = c(0:nrow(clust_profiles)),
                                       k = k_sub, replace = TRUE)
  }

  # Remove all those combinations that only include 1 or 2 cells
  comb <- comb[rowSums(comb != 0) > 2, ]

  # Create all possible combinations
  ## Initialize matrix for increased speed so that it doesn't need to create indexes on the fly
  tmp_mtrx <- matrix(nrow = nrow(comb), ncol = ncol(clust_profiles))
  tmp_metadata <- matrix(nrow = nrow(comb), ncol = nrow(clust_profiles))
  colnames(tmp_metadata) <- rownames(clust_profiles)

  if (verbose) print("Creating synthetic spots"); st_syn_spot <- Sys.time()
  if (verbose) pb_for <- txtProgressBar(min = 0, max = nrow(comb), style = 3) # Progress bar

  for (i in seq_len(nrow(comb))) {
    # Get how many cells of each type we have
    tt <- table(comb[i, ][comb[i, ] != 0])
    tmp_metadata[i, as.numeric(names(tt))] <- tt

    # Add all the profiles together
    row_i <- lapply(names(tt), function(nm) {
      tmp_vec <- tt[[nm]] * clust_profiles[rownames(clust_profiles)[[as.numeric(nm)]], ]
      return(tmp_vec)
    }) %>% purrr::reduce(., `+`)

    # Save mean of the profiles
    tmp_mtrx[i, ] <- row_i / sum(tt)
    # update progress bar
    if (verbose) setTxtProgressBar(pb_for, i)
  }; rm(i, tt, row_i)

  tmp_metadata[is.na(tmp_metadata)] <- 0

  if (verbose) close(pb_for)
  if (verbose) print(sprintf("Creation of %s synthetic spot profiles took: %s minutes",
                            nrow(comb),
                            round(difftime(time1 = Sys.time(),
                                           time2 = st_syn_spot,
                                           units = "mins"), 2)))

  return(list(tmp_mtrx, tmp_metadata))
}
