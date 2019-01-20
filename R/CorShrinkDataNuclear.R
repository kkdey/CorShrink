#' @title Nuclear norm penalization of the covariance matrix
#' @description Performs a nuclear norm penalization of the coavriance matrix.
#'
#' @param data_with_missing The samples by features data matrix. May contain NA values.
#' @param alpha The tuning parameter for the gradient descent(GD) iteration update.
#' @param stepsize The stepsize for the gradient descent algorithms.
#' @param max_iter The maximum number of iterations for the GD.
#' @param tol The tolerance level when to stop the iterations.
#' @param verbose If TRUE, the function prints the objective value on each run, which can be used to check
#'                if the objective is decreasing over iterations (as it should be) or not.
#'
#' @examples
#' data("sample_by_feature_data")
#' out = CorShrinkDataNuclear(sample_by_feature_data, stepsize = 1, max_iter = 100)
#' plot(svd(out$estS)$d)
#' corrplot::corrplot(as.matrix(out$estR), diag = FALSE,
#'         col = colorRampPalette(c("blue", "white", "red"))(200),
#'         tl.pos = "td", tl.cex = 0.4, tl.col = "black",
#'         rect.col = "white",na.label.col = "white",
#'         method = "color", type = "upper")
#'
#' @keywords box shrinkage, correlation
#' @import CVXR
#' @importFrom corrplot corrplot
#' @importFrom stats cor sd cov2cor
#' @importFrom Matrix nearPD
#' @export

CorShrinkDataNuclear <- function(data_with_missing,
                                 alpha = 0.01,
                                 stepsize = 1,
                                 max_iter = 1000,
                                 tol = 1e-04,
                                 verbose = TRUE){

  ##################  Building matrix of common samples for pairwise comparisons  ####################

  binary_indicator = matrix(1, nrow(data_with_missing), ncol(data_with_missing))
  binary_indicator[is.na(data_with_missing)]= 0
  common_samples = t(binary_indicator)%*%binary_indicator
  diag(common_samples) = 0

  ####################  Computing pairwise covariance (correlation) matrices  ###########################

  pairwise_cov = cov(data_with_missing, use = "pairwise.complete.obs")
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
  ########################  Optimization (ADAM or Gradient Descent)  ###########################
  Sigma = diag(sigma_vals) %*% diag(dim(data_with_missing)[2]) %*% diag(sigma_vals)
  estimate_cov = diag(sigma_vals) %*% estimate_cor %*% diag(sigma_vals)
  sd_cov = diag(sigma_vals) %*% sd_cor %*% diag(sigma_vals)
  inverse_sd_cov = 1/sd_cov
  diag(inverse_sd_cov) = 0

  kappa <- 1/max(2*abs(inverse_sd_cov)^2)

  for(iter in 1:max_iter){
    objective = sum((Sigma - estimate_cov)^2*inverse_sd_cov^2) + alpha*sum(diag(Sigma))
    if(verbose){
      cat("value of the objective :", objective, "\n")
    }
    gradient = 2*(Sigma - estimate_cov)*inverse_sd_cov^2 + alpha*diag(dim(Sigma)[1])
    Sigma_tilde = Sigma - stepsize*kappa*gradient
    Sigma_tilde = (Sigma_tilde + t(Sigma_tilde))/2
    eig_decomp = eigen(Sigma_tilde)
    Sigma_tilde_tilde = eig_decomp$vectors %*% diag(pmax(eig_decomp$values, 1e-08)) %*% t(eig_decomp$vectors)
    rel_measure = sum((Sigma_tilde_tilde - Sigma)^2)
    Sigma= Sigma_tilde_tilde
    if(rel_measure < tol){
      break
    }
  }
  sigma2 = mean(diag(estimate_cov)- diag(Sigma))
  if(sigma2 > 0){
    Sigma = Sigma + diag(sigma2, nrow(Sigma))
  }
  R = cov2cor(Sigma)
  ll = list("estS" = Sigma, "estR" = R,  "obj" = objective, "noise_var" = sigma2)
  return(ll)
}
