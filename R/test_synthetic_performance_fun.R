#' If you wanna test the performance of your model on synthetic generated test spots you can use this function to benchmark and get a sense of the model's performance.
#'
#' @param test_spots_metadata Object of class matrix containing the ground truth composition of each spot as obtained from the function syn_spot_comb_topic_fun.R, 2nd element.
#' @param spot_composition_mtrx Object of class matrix with the predicted topic probability distributions for each spot.
#' @return This function returns a list with TP, TN, FP, FN and the Jensen-Shannon Divergence index.
#' @export
#' @examples
#'

test_synthetic_performance <- function(test_spots_metadata_mtrx,
                                       spot_composition_mtrx) {
  # Check variables
  if (!is.matrix(test_spots_metadata_mtrx)) stop("ERROR: test_spots_metadata_mtrx must be a matrix object!")
  if (!is.matrix(spot_composition_mtrx)) stop("ERROR: spot_composition_mtrx must be a matrix object!")

  colnames(spot_composition_mtrx) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                x = colnames(spot_composition_mtrx),
                                perl = TRUE)
  colnames(test_spots_metadata_mtrx) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                   x = colnames(test_spots_metadata_mtrx),
                                   perl = TRUE)
  #load required packages
  suppressMessages(require(philentropy))

  ##### Get TRUE JSD between real-predicted proportions #####
  true_jsd_mtrx <- matrix(nrow = nrow(test_spots_metadata_mtrx), ncol = 1)
  tp <- 0; tn <- 0; fp <- 0; fn <- 0
  for (i in seq_len(nrow(test_spots_metadata_mtrx))) {

    # Create matrix to feed to JSD
    x <- rbind(test_spots_metadata_mtrx[i, ],
               spot_composition_mtrx[i, ])

    # Calculate JSD and save it in true_JSD_mtrx
    if(sum(spot_composition_mtrx[i, ]) > 0) {
      true_jsd_mtrx[i, 1] <- suppressMessages(JSD(x = x, unit = "log2",
                                                  est.prob = "empirical"))
    } else {
      true_jsd_mtrx[i, 1] <- 1
    }

    #### Calculate TP-TN-FP-FN ####
    for (index in colnames(test_spots_metadata_mtrx)) {
      if (x[1, index] > 0 & x[2, index] > 0) {
        tp <- tp + 1
      } else if (x[1, index] == 0 & x[2, index] == 0) {
        tn <- tn + 1
      } else if (x[1, index] > 0 & x[2, index] == 0) {
        fn <- fn + 1
      } else if (x[1, index] == 0 & x[2, index] > 0) {
        fp <- fp + 1
      }
    }; rm(index)

  }; rm(i)

  #### Performance metrics ####
  accuracy <- round((tp + tn) / (tp + tn + fp + fn), 2)
  sensitivity <- round(tp / (tp + fn), 2)
  specificity <- round(tn / (tn + fp), 2)
  precision <- round(tp / (tp + fp), 2)
  recall <- round(tp / (tp + fn), 2)
  F1 <- round(2 * ((precision * recall) / (precision + recall)), 2)

  quants_jsd <- round(quantile(matrixStats::rowMins(true_jsd_mtrx,
                                                    na.rm = TRUE),
                               c(0.25, 0.5, 0.75)), 5)

  cat(sprintf("The following summary statistics are obtained:
              Accuracy: %s,
              Sensitivity: %s,
              Specificity: %s,
              precision: %s,
              recall: %s,
              F1 score: %s,
              JSD quantiles: %s[%s-%s]",
              accuracy, sensitivity, specificity, precision, recall, F1,
              quants_jsd[[2]], quants_jsd[[1]], quants_jsd[[3]]), sep = "\n")

  cat("raw statistics are returned in the list - TP, TN, FP, FN, JSD quantiles",
      sep = "\n")
  return(list(TP = tp, TN = tn, FP = fp, FN = fn, JSD = quants_jsd))
}
