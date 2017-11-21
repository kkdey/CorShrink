#' @title Adaptive shrinkage of correlations from a data matrix
#' @description Performs adaptive shrinkage of the sample correlations starting
#' from a data matrix (possibly containing NAs).
#'
#' @param data The samples by features data matrix. May contain NA values.
#' @param sd_boot A Boolean variable indicating if the standard errors of the
#'                Fisher z-scores should be computed via Bootstrap methods or
#'                through asymptotic formulation of the problem.
#' @param thresh_up Upper threshold for correlations. Defaults to 0.99
#' @param thresh_down Lower threshold for correlations. Defaults to -0.99.
#' @param image if TRUE, plots an image of the shrunk and non-shrunk correlation matrices
#' @param tol The tolerance chosen to check how far apart the CorShrink matrix is from the nearest
#'            positive definite matrix before applying PD completion.
#' @param optmethod The optimization method for EM algorithm - can be one of
#'                  two techniques \code{mixEM} (mixture EM) and \code{mixVBEM}
#'                  (mixture Variational Bayes EM) approaches.The default approach
#'                  is \code{mixEM}.
#' @param ash.control The control parameters for adaptive shrinkage
#'
#' @return Returns an adaptively shrunk version of the sample correlations matrix.
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216; doi: http://dx.doi.org/10.1101/038216
#'
#' @keywords adaptive shrinkage, correlation
#' @importFrom reshape2 melt dcast
#' @import Matrix
#' @import ashr
#' @export


CorShrinkData <- function(data,  sd_boot = FALSE,
                          thresh_up = 0.99, thresh_down = - 0.99,
                          image=FALSE, tol=1e-06,
                          optmethod = "mixEM",
                          ash.control = list()){

  cormat <- cor(data, use = "pairwise.complete.obs")

  if (!sd_boot){
    data2 <- data
    data2[which(!is.na(data2))] <- 1
    data2[which(is.na(data2))] <- 0
    nsamp <- t(data2) %*% data2
    nsamp[nsamp <=2] = 0

    if(dim(nsamp)[1] != dim(cormat)[1] | dim(nsamp)[2] != dim(cormat)[2]){
      stop("dimensions of the matrix of complete samples per pair of variables
           does not match with the correlation matrix")
    }

    out <- CorShrinkMatrix(cormat, nsamp, zscore_sd = NULL,
                           thresh_up = thresh_up, thresh_down = thresh_down,
                           image=image, tol=tol,
                           optmethod = optmethod,
                           ash.control = ash.control)
  }else{
    zscore_sd <- bootcorSE_calc(data)
    out <- CorShrinkMatrix(cormat, nsamp=NULL, zscore_sd = zscore_sd,
                           thresh_up = thresh_up, thresh_down = thresh_down,
                           image=image, tol=tol,
                           optmethod = optmethod,
                           ash.control = ash.control)
  }

  return(out)
}
