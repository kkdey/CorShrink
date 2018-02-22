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
#' @param image character. options for plotting the original or the corshrink matrix.
#'              If \code{image = "both"}, then the function outputs both the plot
#'              for original and shrunk correlationmatrix. If \code{image = "original"},
#'              then the function outputs the correlation plot for the original matrix only.
#'              If \code{image = "corshrink"}, then the function outputs the correlation plot
#'              for the CorShrink matrix only.If \code{image = "output"}, then the function
#'              outputs the saved ggplot2 figure without displaying it. Defaults to "both".
#' @param tol The tolerance chosen to check how far apart the CorShrink matrix is from the nearest
#'            positive definite matrix before applying PD completion.
#' @param image.control Control parameters for the image when \code{image = TRUE}.
#'
#' @param report_model  if TRUE, outputs the full adaptive shrinkage output, else outputs the shrunken vector.
#'                      Defaults to FALSE.
#' @param ash.control The control parameters for adaptive shrinkage
#'
#' @return Returns an adaptively shrunk version of the sample correlations matrix.
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216; doi: http://dx.doi.org/10.1101/038216
#'
#' @examples
#' data("sample_by_feature_data")
#' out <- CorShrinkData(sample_by_feature_data, image = "both")
#'
#' @keywords adaptive shrinkage, correlation
#' @importFrom stats cor sd
#' @export
CorShrinkData <- function(data,  sd_boot = FALSE,
                          thresh_up = 0.99, thresh_down = - 0.99,
                          image = c("both", "original", "corshrink", "output"),
                          tol=1e-06,
                          image.control = list(),
                          report_model = FALSE,
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
                           image = image,
                           tol=tol,
                           image.control = image.control,
                           report_model = report_model,
                           ash.control = ash.control)
  }else{
    zscore_sd <- bootcorSE_calc(data)
    out <- CorShrinkMatrix(cormat, nsamp=NULL, zscore_sd = zscore_sd,
                           thresh_up = thresh_up, thresh_down = thresh_down,
                           image = image,
                           tol=tol,
                           image.control = image.control,
                           report_model = report_model,
                           ash.control = ash.control)
  }

  return(out)
}
