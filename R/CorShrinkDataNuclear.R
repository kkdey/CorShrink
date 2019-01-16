#' @title Nuclear norm penalization of the covariance matrix
#' @description Performs a nuclear norm penalization of the coavriance matrix.
#'
#' @param data_with_missing The samples by features data matrix. May contain NA values.
#' @param method Either gradient descent (GD) or ADAM method of optimization.
#' @param alpha The tuning parameter for the gradientdescent/ADAM iteration update.
#' @param stepsize The stepsize for the ADAM or gradient descent algorithms.
#' @param max_iter The maximum number of iterations for the ADAM/GD.
#' @param tol The tolerance level when to stop the iterations.
#' @param verbose If TRUE, the function prints the objective value on each run, which can be used to check
#'                if the objective is decreasing over iterations (as it should be) or not.
#'
#' @examples
#' data("sample_by_feature_data")
#' out = CorShrinkDataNuclear(sample_by_feature_data, method = "GD", stepsize = 1e-04, max_iter = 50)
#' corrplot::corrplot(as.matrix(out$estR), diag = FALSE,
#'         col = colorRampPalette(c("blue", "white", "red"))(200),
#'         tl.pos = "td", tl.cex = 0.4, tl.col = "black",
#'         rect.col = "white",na.label.col = "white",
#'         method = "color", type = "upper")
#' out = CorShrinkDataNuclear(sample_by_feature_data, method = "ADAM", stepsize = 1e-2, max_iter = 50)
#' corrplot::corrplot(as.matrix(out$estR), diag = FALSE,
#'         col = colorRampPalette(c("blue", "white", "red"))(200),
#'         tl.pos = "td", tl.cex = 0.4, tl.col = "black",
#'         rect.col = "white",na.label.col = "white",
#'         method = "color", type = "upper")
#' @keywords box shrinkage, correlation
#' @import CVXR
#' @importFrom corrplot corrplot
#' @importFrom stats cor sd cov2cor
#' @importFrom Matrix nearPD
#' @export

CorShrinkDataNuclear <- function(data_with_missing,
                                 method,
                                 alpha = 0.01,
                                 stepsize = 1e-03,
                                 max_iter = 1000,
                                 tol = 1e-04,
                                 verbose = TRUE){

  #########################  Defining the method as default   ##############################

  if(missing(method)){
    method = "GD"
    message("By default using the gradient descent method of optimization, else put method = 'ADAM' for accelerated gradient descent")
  }
  if(method != "ADAM" && method != "GD"){
    stop("The method must be either ADAM or GD (gradient descent)")
  }
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

  mt= 0
  vt = 0
  beta1=0.9
  beta2 = 0.99

  for(iter in 1:max_iter){
    objective = sum((Sigma - estimate_cov)^2*inverse_sd_cov^2) + alpha*sum(diag(Sigma))
    if(verbose){
      cat("value of the objective :", objective, "\n")
    }
    gradient = 2*(Sigma - estimate_cov)*inverse_sd_cov^2 + alpha*diag(dim(Sigma)[1])
    if(method == "ADAM"){
      mt = beta1*mt+ (1-beta1)*gradient
      vt = beta2*vt + (1-beta2)*gradient^2
      mthat = mt/(1-beta1)
      vthat = vt/(1-beta2)
      Sigma_tilde = Sigma - (stepsize/(sqrt(vthat)+1e-08))*mthat
      Sigma_tilde_tilde = as.matrix(Matrix::nearPD(Sigma_tilde)$mat)
    }else if(method == "GD"){
      Sigma_tilde = Sigma - stepsize*gradient
      Sigma_tilde_tilde = as.matrix(Matrix::nearPD(Sigma_tilde)$mat)
    }
    rel_measure = sum((Sigma_tilde_tilde - Sigma)^2)
    Sigma= Sigma_tilde_tilde
    if(rel_measure < tol){
      break
    }
  }
  R = cov2cor(Sigma)
  ll = list("estS" = Sigma, "estR" = R,  "obj" = objective)
  return(ll)
}
