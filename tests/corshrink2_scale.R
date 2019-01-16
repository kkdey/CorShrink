

#######################  CorShrink2Scale with nuclear norm  #################################

library(CorShrink)
data("sample_by_feature_data")
data("common_samples")

#####################   Compute pairwise covariance matrix  ############################

pairwise_cov = cov(sample_by_feature_data, use = "pairwise.complete.obs")
sigma_vals = sqrt(diag(pairwise_cov))
pairwise_cor = cov2cor(pairwise_cov)
pairwise_cor[is.na(pairwise_cor)] = 0
pairwise_cor[pairwise_cor > 0.95] = 0.95
pairwise_cor[pairwise_cor < -0.95] = -0.95
diag(pairwise_cor) = 1
common_samples[common_samples <= 2] = 2
pairwise_cov = diag(sigma_vals) %*% pairwise_cor %*% diag(sigma_vals)

#####################  Sample size bias factor   ###########################

samples_adjustment = 1/(common_samples - 1)  + 2/(common_samples - 1)^2

####################  Bias and scalings   ####################################

bias_cor = samples_adjustment*pairwise_cor*(1 + pairwise_cor)*(1 - pairwise_cor)
sd_cor = sqrt(samples_adjustment)*(1+pairwise_cor)*(1-pairwise_cor)

estimate_cor = pairwise_cor + bias_cor
inverse_sd_cor = 1/sd_cor
diag(inverse_sd_cor) = 0

R = diag(dim(sample_by_feature_data)[2])
max_iter = 10000
lambda = 0.0001
step = 1
mt= 0
vt = 0
beta1=0.9
beta2 = 0.99
method = "ADAM"
for(iter in 1:max_iter){
  #objective = sum(((R - estimate_cor)*inverse_sd_cor)^2)
  objective = sum((R - estimate_cor)^2*inverse_sd_cor^2) + lambda*sum(diag(R))
  cat("value of the objective :", objective, "\n")
  #gradient  = 2*(R - estimate_cor)*(inverse_sd_cor^2)
  gradient = 2*(R - estimate_cor)*inverse_sd_cor^2 + lambda*diag(dim(R)[1])
  if(method == "ADAM"){
    mt = beta1*mt+ (1-beta1)*gradient
    vt = beta2*vt + (1-beta2)*gradient^2
    mthat = mt/(1-beta1)
    vthat = vt/(1-beta2)
    R_tilde = R - (stepsize/(sqrt(vthat)+1e-08))*mthat
    R_tilde_tilde = cov2cor(as.matrix(nearPD(R_tilde)$mat))
  }else if(method == "GD"){
    R_tilde = R - stepsize*gradient
    R_tilde_tilde = cov2cor(as.matrix(nearPD(R_tilde)$mat))
  }
  col2 <- c("blue", "white", "red")
  rel_measure = sum((R_tilde_tilde - R)^2)
  cat("the relative error in approximation between two steps: ", rel_measure, "\n")
  R= R_tilde_tilde
}


col2 <- c("blue", "white", "red")
corrplot(as.matrix(R_tilde_tilde), diag = FALSE, col = colorRampPalette(col2)(200), tl.pos = "td",
         tl.col = "black", tl.cex = 0.8, rect.col = "white",
         na.label.col = "white", method = "color", type = "upper")





R <- Semidef(dim(pairwise_cor)[1])
obj <- Minimize(p_norm(R, norm1) + tau*p_norm(mul_elemwise(inverse_sd_cor, square(R - estimate_cor)), 1))
constraints = list(diag(R) == 1)
prob <- Problem(obj, constraints)
result <- solve(prob)
R_hat = as.matrix(result$getValue(R))


