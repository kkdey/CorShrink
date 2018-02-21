

#############  test corshrink_word2vec   #############################

model_list <- list()

for(ll in 1:50){
  model_list[[ll]] <- read.vectors(paste0("../internal_data/fake_nation/", ll, ".pool.bin"))
}
