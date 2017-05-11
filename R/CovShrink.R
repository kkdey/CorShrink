#' @title Adaptive shrinkage of a covariance matrix through shrinking of correlations.
#'
#' @description Performs adaptive shrinkage of covariance matrix through shrinkage of the
#' underlying correlation matrix using wither of the three schemes \code{CorShrinkML},
#' \code{CorShrinkVEM} or \code{CorShrinkVEM2}.
#'
#' @param data The samples by features data matrix.
#' @param nsamples The number of samples used for estimating correlation matrix.
#' @param nullweight The weight on the null component for the mixture normal prior used. Defaults to 10.
#' @param null.comp The number of null components chosen. Defaults to 1.
#' @param type The type of adaptive correlation matrix shrinkage to be used.
#' @param ash.control The control parameters for adaptive shrinkage
#'
#' @return Returns a shrunk version of the sample covariance matrix (the output is also a covariance matrix)
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216; doi: http://dx.doi.org/10.1101/038216
#' @keywords adaptive shrinkage, correlation
#' @importFrom reshape2 melt dcast
#' @import Matrix
#' @import ashr
#' @export


CovShrink <- function(data, nsamples, sd_boot = FALSE, partial = FALSE,
                      type = c("ML", "VEM"), nboot = 50,
                      nullweight = 10, null.comp = 1, thresh_up = 0.999,
                      thresh_down = 0.001, tol = 1e-06, ash.control = list()){
  cov_sample <- cov(data);
  cor_sample <- cov2cor(cov_sample);
  if(partial == TRUE){
    cor_sample <- corpcor::cor2pcor(cor_sample)
  }
  if(sd_boot){
    cor_transform_sd_vec <- bootcorSE_calc(data, nboot=nboot)
  }else{
    cor_transform_sd_vec <- NULL
  }
  if(type=="ML"){
    cor.shrunk <- CorShrinkML(cor_sample, nsamples, sd_boot = sd_boot,
                              cor_transform_sd_vec = cor_transform_sd_vec,
                              image = FALSE, thresh_up = thresh_up,
                              thresh_down = thresh_down, tol=tol,
                              ash.control = ash.control)$ash_cor_PD
  }else if (type == "VEM"){
    cor.shrunk <- CorShrinkVEM(cor_sample, nsamples,
                               sd_boot = sd_boot,
                               cor_transform_sd_vec = cor_transform_sd_vec,
                               nullweight = nullweight,
                               null.comp =null.comp,
                               thresh_up = thresh_up,
                               thresh_down = thresh_down, tol=tol)$ash_cor_PD
  }else{
    stop("the type provided must be one of ML and VEM")
  }

  diags <- diag(cov_sample)
  diags[diags < 0] = 0
  cov.shrunk <- diag(sqrt(diags))%*%cor.shrunk%*%diag(sqrt(diags))
  return(cov.shrunk)
}

