
##############################  test  1   ###################################


#######################  test CorShrinkMatrix   #######################

common_samples_mat <- get(load(file = "../data/common_samples.rda"))
data <- get(load(file = "../data/sample_by_feature_data.rda"))
cormat <- cor(data, use = "pairwise.complete.obs")

data2 <- data
data2[which(!is.na(data2))] <- 1
data2[which(is.na(data2))] <- 0
nsamp <- t(data2) %*% data2
nsamp[nsamp <=2] = 0

if(dim(nsamp)[1] != dim(cormat)[1] | dim(nsamp)[2] != dim(cormat)[2]){
  stop("dimensions of the matrix of complete samples per pair of variables
       does not match with the correlation matrix")
}

out <- CorShrinkMatrix(cormat, nsamp, image = TRUE, optmethod = "mixEM")
out <- CorShrinkMatrix(cormat, nsamp, image = TRUE, optmethod = "mixVBEM")

##################  test CorShrinkData   ###############################


data <- get(load(file = "../data/sample_by_feature_data.rda"))
out <- CorShrinkData(data, image = TRUE, optmethod = "mixEM")


##################   test  CorShrinkVector  ###########################

cor_vec <- c(-0.56, -0.4, 0.02, 0.2, 0.9, 0.8, 0.3, 0.1, 0.4)
nsamp_vec <- c(10, 20, 30, 4, 50, 20, 20, 10, 3)
out <- CorShrinkVector(corvec = cor_vec, nsamp_vec = nsamp_vec,
                       optmethod = "mixVBEM")

########  Previous formulations used the asymptotic distribution of the
########  Fisher Z-scores. Now we use resampling for more finite sample
########  distributions under non-normality assumptions of Fisher z-scores.


##########  Bootstrap standard error of z scores calculator   #############

data <- get(load(file = "../data/sample_by_feature_data.rda"))
zscoreSDmat <- bootcorSE_calc(data)


##########  use Bootstrap standard errors in Corshrink  ##################

out <- CorShrinkData(data, sd_boot = TRUE)
out <- CorShrinkMatrix(cormat, zscore_sd = zscoreSDmat)

####################  Output assessement   #######################

dim(out$ash_cor_only)
dim(out$ash_cor_PD)

out$ash_cor_only[1:5,1:5]
out$ash_cor_PD[1:5,1:5]




