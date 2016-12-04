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
#'
#' @return Returns a shrunk version of the sample covariance matrix (the output is also a covariance matrix)
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216; doi: http://dx.doi.org/10.1101/038216
#' @keywords adaptive shrinkage, correlation
#' @importFrom reshape2 melt dcast
#' @import Matrix
#' @import ashr
#' @export


CovShrink <- function(data, nsamples, nullweight=10, null.comp=1,
                      type = c("ML", "VEM", "VEM2")){
  cov_sample <- cov(data);
  cor_sample <- cov2cor(cov_sample);
  if(type=="ML"){
    cor.shrunk <- CorShrinkML(cor_sample, nsamples)$ash_cor_PD
  }else if (type == "VEM"){
    cor.shrunk <- CorShrinkVEM(cor_sample, nsamples, nullweight = nullweight, null.comp =null.comp)$ash_cor_PD
  }else if (type == "VEM2"){
    cor.shrunk <- CorShrinkVEM2(cor_sample, nsamples, nullweight = nullweight, null.comp =null.comp)$ash_cor_PD
  }else{
    stop("the type provided must be one of ML, VEM and VEM2")
  }

  diags <- eigen(cov_sample)$values
  cov.shrunk <- diag(sqrt(diags))%*%cor.shrunk%*%diag(sqrt(diags))
  return(cov.shrunk)
}

