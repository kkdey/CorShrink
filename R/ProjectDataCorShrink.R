#' @title Project data on the shrunk correlation or the covariance model space.
#'
#' @description Performs projection of the sample data into the the shrunk
#' correlation or the covariance model space, where the shrinkage is performed adaptively
#' on the covariance model using the \code{CovShrink} function which in turn uses one of the
#' three correlation shrinkage methods - \code{CorShrinkML},
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
#' @importFrom expm sqrtm
#' @export
#'

ProjectDataCorShrink <- function(data, nsamples, nullweight=10, null.comp=1,
                                 type = c("ML", "VEM", "VEM2")){
  cov_sample <- cov(data);
  cov_shrunk <- CovShrink(data, nsamples, nullweight = nullweight,
                          null.comp = null.comp, type = type)
  projected_dat <- expm::sqrtm(cov_shrunk)%*% expm::sqrtm(solve(cov_sample)) %*% data
  return(projected_dat)
}
