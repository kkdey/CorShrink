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
#' @param image_original if TRUE, plots an image of the non-shrunk original matrix
#'                       of correlations.
#' @param image_corshrink if TRUE, plots an image of the shrunk matrix
#'                       of correlations.
#' @param tol The tolerance chosen to check how far apart the CorShrink matrix is from the nearest
#'            positive definite matrix before applying PD completion.
#' @param image.control Control parameters for the image when \code{image = TRUE}.
#' @param optmethod The optimization method for EM algorithm - can be one of
#'                  two techniques \code{mixEM} (mixture EM) and \code{mixVBEM}
#'                  (mixture Variational Bayes EM) approaches.The default approach
#'                  is \code{mixEM}.
#' @param report_model  if TRUE, outputs the full adaptive shrinkage output, else outputs the shrunken vector.
#'                      Defaults to FALSE.
#' @param ash.control The control parameters for adaptive shrinkage
#'
#' @return Returns an adaptively shrunk version of the sample correlations matrix.
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216; doi: http://dx.doi.org/10.1101/038216
#'
#' @examples
#' data <- get(load(system.file("extdata", "sample_by_feature_data.rda",
#'                    package = "CorShrink")))
#' out <- CorShrinkData(data, sd_boot = FALSE, image_original = TRUE,
#'                       image_corshrink = TRUE, optmethod = "mixEM",
#'                       image.control = list(x.cex = 0.3, y.cex = 0.3))
#'
#' @keywords adaptive shrinkage, correlation
#' @importFrom reshape2 melt dcast
#' @importFrom stats cor sd
#' @importFrom utils modifyList
#' @export
CorShrinkData <- function(data,  sd_boot = FALSE,
                          thresh_up = 0.99, thresh_down = - 0.99,
                          image_original=FALSE, image_corshrink = FALSE,
                          tol=1e-06,
                          image.control = list(),
                          optmethod = "mixEM",
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
                           image_original=image_original,
                           image_corshrink = image_corshrink,
                           tol=tol,
                           image.control = image.control,
                           optmethod = optmethod,
                           report_model = report_model,
                           ash.control = ash.control)
  }else{
    zscore_sd <- bootcorSE_calc(data)
    out <- CorShrinkMatrix(cormat, nsamp=NULL, zscore_sd = zscore_sd,
                           thresh_up = thresh_up, thresh_down = thresh_down,
                           image_original=image_original,
                           image_corshrink = image_corshrink,
                           tol=tol,
                           image.control = image.control,
                           optmethod = optmethod,
                           report_model = report_model,
                           ash.control = ash.control)
  }

  return(out)
}
