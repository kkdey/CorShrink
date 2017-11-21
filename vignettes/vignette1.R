
##############################  test  1   ###################################


common_samples_mat <- get(load(file = "../data/common_samples.rda"))
data <- get(load(file = "../data/sample_by_feature_data.rda"))
cormat <- pcor(data, use = "pairwise.complete.obs")

data2 <- data
data2[which(!is.na(data2))] <- 1
data2[which(is.na(data2))] <- 0
nsamp <- t(data2) %*% data2
nsamp[nsamp <=2] = 0

if(dim(nsamp)[1] != dim(cormat)[1] | dim(nsamp)[2] != dim(cormat)[2]){
  stop("dimensions of the matrix of complete samples per pair of variables
       does not match with the correlation matrix")
}



out <- CorShrinkMatrix(cormat, nsamp, image = TRUE)

