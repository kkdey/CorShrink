#' @title CorShrink on word rankings in word2vec model
#'
#' @description This function computes word rankings on the word2vec model
#' from CorShrink using resampling methods.The method is similar to the adaptive shrinkage method for modeling false discovery rates proposed
#' in Stephens 2016 (see reference).
#'
#'
#' @param model_true the true word2vec model
#' @param model_boot_list a list of bootstrap word2vec models with resampling lines of text.
#' @param word_vec A vector of words whose rankings are to be modified
#' @param num_related_words Number of related words to consider per word in \code{word_vec}
#' @param ash.control The control parameters for adaptive shrinkage
#'
#' @return Returns adaptively shrunk cosine similarity values of the words of interest in \code{word_vec}.
#'
#' @keywords adaptive shrinkage, correlation
#' @import ashr
#' @export

CorShrink_word2vec <- function(model_true, model_boot_list,
                               word_vec, num_related_words = 50,
                               ash.control = list()){

  ash.control.default = list(pointmass = TRUE, prior = "nullbiased", gridmult = 2,
                             mixcompdist = "halfuniform", nullweight = 10,
                             outputlevel = 2, fixg = FALSE, optmethod="mixEM")
  ash.control <- modifyList(ash.control.default, ash.control)

  ############  We initialize the list which we will return at end containing the results
  ##########  from the analysis

   return_list <- vector(mode = "list", length(word_vec))

  ###########  the number of bootstrap samples : equal to length of model_boot_list  #####


    numboot <- length(model_boot_list)

  ############  for each word, carry out Corshrink over related words  ##################

    for(numvec in 1:length(word_vec)){
      ############  taking a word from the word vector passed
         word <- word_vec[numvec]

         ############  finding the neighbors or related words to that word  #############


         word_nbrs <- as.character()
         for(ll in 1:numboot){
           word_nbrs <- c(word_nbrs, names(nearest_to(model_boot_list[[ll]],
                                                      model_boot_list[[ll]][[paste0(word)]], num_related_words)))
         }

         listed_words <- unique(word_nbrs)
         listed_words <- listed_words[-1]

         ############ sd of cosine-similarity of related words to the given word
         ############  from bootstrap samples   #################

         cosine_transform_sd_vec <- array(0, length(listed_words))
         cosine_boot_extended <- array(0, length(listed_words))

         mat_temp <- matrix(0, numboot, length(listed_words))

         for(i in 1:numboot){
           cosine_boot <- unlist(as.data.frame(cosineSimilarity(model_boot_list[[i]][[word, average = FALSE]],model_boot_list[[i]][[listed_words, average = FALSE]])))
           match_indices <- match(names(cosine_boot), listed_words)
           cosine_boot_extended[match_indices] <- cosine_boot
           cosine_transform_boot_extended <- 0.5*log((1+cosine_boot_extended)/(1-cosine_boot_extended))
           mat_temp[i,] <- cosine_transform_boot_extended
        }

         cosine_transform_sd_vec <- apply(mat_temp, 2, function(x) return(sd(x)))


         cosine_est_extended <- array(0, length(listed_words))

         cosine_est <- unlist(as.data.frame(cosineSimilarity(model_true[[word, average = FALSE]], model_true[[listed_words, average = FALSE]])))
         match_indices <- match(names(cosine_est), listed_words)
         cosine_est_extended[match_indices] <- cosine_est
         cosine_transform_est_extended <- 0.5*log((1+cosine_est_extended)/(1-cosine_est_extended))
         names(cosine_est_extended) <- listed_words

         ash_out <- do.call(ashr::ash, append(list(as.vector(cosine_transform_est_extended),
                              as.vector(cosine_transform_sd_vec)),
                              ash.control))
         fitted_cos <- (exp(2*ash_out$result$PosteriorMean)-1)/(exp(2*ash_out$result$PosteriorMean)+1)
         names(fitted_cos) <- listed_words
         return_list[[numvec]] = list(similar_words = listed_words,
                                      cosine_est = cosine_est_extended,
                                      sd_cosine_transform_est = cosine_transform_sd_vec,
                                      ash_result = ash_out,
                                      ash_cosine_est = fitted_cos)
         cat("Processed the ash shrinkage of word2vec for word: ", paste0(word), "\n")
    }
    return(return_list)
}
