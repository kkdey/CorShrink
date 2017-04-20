
######################  CorShrink.word2vecc   ##########################################
##########   corshrink on word2vec models fitted by the WordVectors package  ###############


CorShrink_word2vec <- function(model_true, model_boot_list,
                               word_vec, num_related_words = 50){

  ############  We initialize the list which we will return at end containing the results
  ##########  from the analysis

   return_list <- vector(mode = "list", length(word_vec))

  ###########  the numbe of bootstrap samples : equal to length of model_boot_list  #####


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
         

         # cosine_transform_sd_list <- parallel::mclapply(1:length(listed_words),
         #                                    function(l){
         #                                             cosine_boot <- array(0, numboot)
         #                                             for(i in 1:numboot){
         #                                                  cosine_boot[i] <- model_boot_list[[i]][[paste0(word)]] %>% cosineSimilarity(model_boot_list[[i]][[paste0(listed_words[l])]])
         #                                           }
         #                                              cosine_boot[which(cosine_boot == "NaN")] = 0
         #                                              cosine_transform_boot <- 0.5*log((1+cosine_boot)/(1-cosine_boot))
         #                                              cosine_transform_sd <- sd(cosine_transform_boot)
         #                                              return(cosine_transform_sd)
         #                                          }, mc.cores = parallel::detectCores())
         # 
         # cosine_transform_sd_vec <- unlist(cosine_transform_sd_list)


#
#          for(m in 1:length(cosine_transform_sd)){
#            cosine_boot <- array(0, numboot)
#            for(i in 1:numboot){
#              cosine_boot[i] <- model_boot_list[[i]][[paste0(word)]] %>% cosineSimilarity(model_boot_list[[i]][[paste0(listed_words[m])]])
#            }
#            cosine_boot[which(cosine_boot == "NaN")] = 0
#            cosine_transform_boot <- 0.5*log((1+cosine_boot)/(1-cosine_boot))
#            cosine_transform_sd[m] <- sd(cosine_transform_boot)
#          }


         ##########  cosine similarity and its transform
         ##########  estimation for the original data  #################

         cosine_est_extended <- array(0, length(listed_words))
         
         cosine_est <- unlist(as.data.frame(cosineSimilarity(model_true[[word, average = FALSE]], model_true[[listed_words, average = FALSE]])))
         match_indices <- match(names(cosine_est), listed_words)
         cosine_est_extended[match_indices] <- cosine_est         
         cosine_transform_est_extended <- 0.5*log((1+cosine_est_extended)/(1-cosine_est_extended))
         names(cosine_est_extended) <- listed_words
         
         # for(m in 1:length(cosine_est)){
         #   cosine_est[m] <- model_true[[paste0(word)]] %>% cosineSimilarity(model_true[[paste0(listed_words[m])]])
         #   cosine_est[which(cosine_est == "NaN")] = 0
         # }


         ##########   Applying the ash shrinkage on cosine transform and inverting back
         #########  the posterior

         ash_out <- ashr::ash(as.vector(cosine_transform_est_extended), as.vector(cosine_transform_sd_vec))
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
