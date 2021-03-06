

########################   Test iCorShrink2Data  on PLIN1 gene correlation matrix  ##########################

library(CorShrink)
data("sample_by_feature_data")
alpha = 0.1
lambda = 0.8
data_with_missing = sample_by_feature_data


library(CVXR)
if(missing(alpha)){
  stop("The value of alpha as in shrinkage intensity not provided, please specify alpha")
}

##################  Building matrix of common samples for pairwise comparisons  ####################

binary_indicator = matrix(1, nrow(data_with_missing), ncol(data_with_missing))
binary_indicator[is.na(data_with_missing)]= 0
common_samples = t(binary_indicator)%*%binary_indicator
diag(common_samples) = 0


###############  compute pairwise covariance/correlation matrix   #######################

pairwise_cov = cov(data_with_missing, use = "pairwise.complete.obs")
sigma_vals = sqrt(diag(pairwise_cov))
pairwise_cor = cov2cor(pairwise_cov)
pairwise_cor[is.na(pairwise_cor)] = 0
pairwise_cor[pairwise_cor > 0.95] = 0.95
pairwise_cor[pairwise_cor < -0.95] = -0.95
diag(pairwise_cor) = 1
common_samples[common_samples <= 2] = 2
pairwise_cov = diag(sigma_vals) %*% pairwise_cor %*% diag(sigma_vals)

###############  Compute pairwise Fisher Z-scores   ##########################

pairwise_zscores = apply(pairwise_cor, c(1,2), function(x) return (0.5*log((1+x)/(1-x))))
diag(pairwise_zscores) = 0

################  Bound on the covariances   ##########################

bound1 = 12*exp(2*pairwise_zscores)/((exp(2*pairwise_zscores) + 1)^2)
zscores_sd_1 = sqrt(1/(common_samples - 1) + 2/(common_samples - 1)^2)
overall_bound_1 = bound1*zscores_sd_1 + zscores_sd_1^2*2*sqrt(3)

common_samples_2 = matrix(dim(data_with_missing)[1], dim(common_samples)[1], dim(common_samples)[2])
zscores_sd_2 = sqrt(1/(common_samples_2 - 1) + 2/(common_samples_2 - 1)^2)
overall_bound_2 = bound1*zscores_sd_2 + zscores_sd_2^2*2*sqrt(3)

delta = apply(abs(overall_bound_1) + abs(overall_bound_2), c(1,2), function(x) return(pmin(2,x)))
diag(delta) = 0
delta_cov = diag(sigma_vals) %*% delta %*% diag(sigma_vals)


#################################  Use adaptive regularized GLASSO  ############################

MAX_ITER = 10
Weight = matrix(1, nrow(delta_cov), ncol(delta_cov))
Omega_hat = diag(nrow(delta_cov))
rel_error_vec = c()
power = 0.05
for(iter in 1:MAX_ITER){
  Omega_hat_old = Omega_hat
  Omega <- Semidef(dim(common_samples)[1])
  scale <- Weight*abs(alpha) + lambda*delta_cov
  obj = Minimize(-log_det(Omega) + matrix_trace(Omega %*% pairwise_cov) + sum(mul_elemwise(scale, abs(Omega))))
  prob <- Problem(obj)
  result <- solve(prob)
  R_hat <- base::solve(result$getValue(Omega))
  Omega_hat <- result$getValue(Omega)
  rel_error = sqrt(sum((Omega_hat - Omega_hat_old)^2))/sqrt(sum(Omega_hat^2))
  cat("The relative error (Frobenius) between last two runs is", rel_error, "\n")
  rel_error_vec = c(rel_error_vec, rel_error)
  Weight = 1/((abs(Omega_hat))^{power} + 0.01)
  if(rel_error < 1e-03){
    cat("Finished iterations!")
    break
  }
}

PR_hat = -cov2cor(as.matrix(Omega_hat))
diag(PR_hat) = 1


