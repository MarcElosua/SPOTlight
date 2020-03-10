#' This function processes the data and trains the LDA model
#'
#' @param lda_mod Object of class LDA_Gibbs.
#' @param k topic we want to visualize.
#' @param n_terms number of terms to show.
#' @return The n_terms most likely genes defining that cluster.
#' @export
#' @examples
#'

topic_viz <- function(lda_mod,
                      k,
                      n_terms) {

  if (is(lda_mod)[[1]] != "LDA_Gibbs") stop("ERROR: lda_mod must be an LDA_Gibbs object!")
  if (!is.numeric(k)) stop("ERROR: k must be an integer!")
  if (!is.numeric(n_terms)) stop("ERROR: n_terms must be an integer!")

  suppressMessages(require(wordcloud2))
  suppressMessages(require(topicmodels))

  tmp_result <- topicmodels::posterior(lda_mod)

  # visualize topics as word cloud
  topic_to_viz <- k # change for your own topic of interest

  # select to 40 most probable terms from the topic by sorting the term-topic-probability vector in decreasing order
  top_n_terms <- sort(tmp_result$terms[topic_to_viz, ],
                     decreasing = TRUE)[1:n_terms]
  words <- names(top_n_terms)

  # extract the probabilites of each of the 40 terms
  probabilities <- sort(tmp_result$terms[topic_to_viz, ],
                        decreasing = TRUE)[1:n_terms]

  # visualize the terms as wordcloud
  wordcloud2::wordcloud2(data.frame(words, probabilities),
                       shuffle = FALSE, size = 0.8)
  return(top_n_terms)
}
