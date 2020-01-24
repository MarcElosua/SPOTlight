#' If you wanna test the performance of your model on synthetic generated test spots you can use this function to benchmark and get a sense of the model's performance.
#'
#' @param test_spots_metadata Object of class matrix with the predicted topic probability distributions for each spot. Output from the lda_prediction function.
#' @param spot_composition_mtrx list obtained from the function syn_spot_comb_topic_fun.R. The 1st element is a matrix with topic profiles of all the synthetic spots generated, the 2nd element is the composition of each synthetic spot.
#' @return This function returns a list with TP, TN, FP, FN and the Jensen-Shannon Divergence index.
#' @export
#' @examples
#' 

test_synthetic_performance <- function(test_spots_metadata_mtrx, spot_composition_mtrx){
  # Check variables
  if( !( is.matrix(test_spots_metadata_mtrx) | is.data.frame(prediction) ) ) {stop("ERROR: prediction must be a matrix object!")}
  if( !( is.list(spot_composition_mtrx) | is.data.frame(prediction) ) ) {stop("ERROR: syn_spots_ls must be the list obtained from the function syn_spot_comb_topic_fun().")}
  
  #load required packages
  suppressMessages(require(philentropy))
  
  ##### Get TRUE JSD between real-predicted proportions #####
  true_JSD_mtrx <- matrix(nrow = nrow(test_spots_metadata_mtrx), ncol = 1)
  tp=0;tn=0;fp=0;fn=0
  for (i in 1:nrow(spot_composition_mtrx)) {
    # Create matrix to feed to JSD
    x <- rbind('truth'=test_spots_metadata_mtrx[i,],'pred'=spot_composition_mtrx[i,])
    
    # Calculate JSD and save it in true_JSD_mtrx
    true_JSD_mtrx[i,1] <- JSD(x = x, unit = "log2", est.prob = 'empirical')  
    
    #### Calculate TP-TN-FP-FN ####
    for (index in c(1:ncol(test_spots_metadata_mtrx))) {
      if(sum(x[1,index]) > 0 & sum(x[2,index]) > 0) {tp=tp+1
      } else if (sum(x[1,index]) == 0 & sum(x[2,index]) == 0) {tn=tn+1
      } else if (sum(x[1,index]) > 0 & sum(x[2,index]) == 0) {fn=fn+1
      } else if (sum(x[1,index]) == 0 & sum(x[2,index]) > 0) fp=fp+1
    }; rm(index)
    
  }; rm(i)
  
  #### Performance metrics ####
  accuracy = round((tp+tn)/(tp+tn+fp+fn),2)
  sensitivity = round(tp/(tp+fp),2)
  specificity = round(tn/(tn+fn),2)
  precision = round(tp/(tp+fp),2)
  recall = round(tp/(tp+fn),2)
  F1 = round(2*((precision*recall)/(precision+recall)),2)
  
  quants_JSD <- round(quantile(matrixStats::rowMins(true_JSD_mtrx,na.rm = TRUE),c(0.25,0.5,0.75)),5)
  
  cat(sprintf('The following summary statistics are obtained:
              Accuracy: %s,
              Sensitivity: %s,
              Specificity: %s,
              precision: %s,
              recall: %s,
              F1 score: %s,
              JSD quantiles: %s[%s-%s]',
              accuracy, sensitivity, specificity, precision, recall, F1, quants_JSD[[2]], quants_JSD[[1]], quants_JSD[[3]]))
  
  print('raw statistics are returned in the list - TP, TN, FP, FN, JSD quantiles')
  return(list(TP=tp,TN=tn,FP=fp,FN=fn,JSD=quants_JSD))
}