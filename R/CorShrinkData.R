#' @title Adaptive shrinkage of correlations from a data matrix
#' @description Performs adaptive shrinkage of the sample correlations starting
#' from a data matrix (possibly containing NAs).
#'
#' @param data The samples by features data matrix. May contain NA values.
#' @param sd_boot A Boolean variable indicating if the standard errors of the
#'                Fisher z-scores should be computed via Bootstrap methods or
#'                through asymptotic formulation of the problem.
#' @param type character. Either "cor" or "pcor" - depending on whether to use
#'             correlation or partial correlation. Default is "cor".
#' @param cor_method The method of correlation used. May be "pearson", "spearman" or "kendall"
#'                   depending on the type of correlation to be used by the user.
#' @param thresh_up Upper threshold for correlations. Defaults to 0.99
#' @param thresh_down Lower threshold for correlations. Defaults to -0.99.
#' @param image character. options for plotting the original or the corshrink matrix.
#'              If \code{image = "both"}, then the function outputs both the plot
#'              for original and shrunk correlationmatrix. If \code{image = "original"},
#'              then the function outputs the correlation plot for the original matrix only.
#'              If \code{image = "corshrink"}, then the function outputs the correlation plot
#'              for the CorShrink matrix only.If \code{image = "output"}, then the function
#'              outputs the saved ggplot figure without displaying it. If \code{image = "null"},
#'              no image is output.Defaults to "both".
#' @param tol The tolerance chosen to check how far apart the CorShrink matrix is from the nearest
#'            positive definite matrix before applying PD completion.
#' @param nboot The number of bootstrap samples if \code{sd_boot = TRUE}.
#' @param image.control Control parameters for the image when \code{image = TRUE}.
#' @param report_model  if TRUE, outputs the full adaptive shrinkage output,
#'                      else outputs the shrunken vector. Defaults to FALSE.
#' @param maxiter The maximum number of iterations run for the adaptive shrinkage EM algorithm.
#'                 Default is 1000.
#' @param ash.control The control parameters for adaptive shrinkage
#'
#' @return Returns an adaptively shrunk version of the sample correlations matrix.
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216; doi: http://dx.doi.org/10.1101/038216
#'
#' @examples
#' data("sample_by_feature_data")
#' out <- CorShrinkData(sample_by_feature_data, image = "null")
#'
#' @keywords adaptive shrinkage, correlation
#' @importFrom stats cor sd
#' @importFrom corpcor cor2pcor
#' @export
CorShrinkData <- function(data,  sd_boot = FALSE,
                          type = "cor",
                          cor_method,
                          thresh_up = 0.99, thresh_down = - 0.99,
                          image = c("both", "original", "corshrink", "output", "null"),
                          tol=1e-06,
                          dosym=TRUE,
                          nboot = 50,
                          image.control = list(),
                          report_model = FALSE,
                          maxiter = 1000,
                          ash.control = list()){

  if(missing(cor_method)){
    cor_method = "pearson"
  }else{
    if(length(cor_method) > 1){
      stop("cor_method must be either `pearson`, `kendall` or `spearman`")
    }else{
      if(!(cor_method %in% c("pearson", "spearman", "kendall"))){
        stop("cor_method must be either `pearson`, `kendall` or `spearman`")
      }
    }
  }

  if(!(type %in% c("cor", "pcor"))){
    stop("type can either be `cor` or `pcor`")
  }

  cormat <- cor(data, use = "pairwise.complete.obs", method = cor_method)

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

    if(type == "cor"){
      out <- CorShrinkMatrix(cormat, nsamp, zscore_sd = NULL,
                             thresh_up = thresh_up, thresh_down = thresh_down,
                             image = image,
                             tol=tol,
                             dosym=dosym,
                             image.control = image.control,
                             report_model = report_model,
                             maxiter = maxiter,
                             ash.control = ash.control)
    }else if (type == "pcor"){
      out <- CorShrinkMatrix(cormat, nsamp, zscore_sd = NULL,
                             thresh_up = thresh_up, thresh_down = thresh_down,
                             image = image,
                             tol=tol,
                             dosym=dosym,
                             image.control = image.control,
                             report_model = report_model,
                             maxiter = maxiter,
                             ash.control = ash.control)
      out$cor <- corpcor::cor2pcor(out$cor)
    }
  }else{
    if(type == "pcor"){
      zscore_sd <- bootpcorSE_calc(data, cor_method = cor_method, nboot = nboot)
      pcormat <- corpcor::cor2pcor(cormat)
      out <- CorShrinkMatrix(pcormat, nsamp=NULL, zscore_sd = zscore_sd,
                             thresh_up = thresh_up, thresh_down = thresh_down,
                             image = image,
                             tol=tol,
                             dosym=dosym,
                             image.control = image.control,
                             report_model = report_model,
                             maxiter = maxiter,
                             ash.control = ash.control)
    }else if (type == "cor"){
      zscore_sd <- bootcorSE_calc(data, cor_method = cor_method, nboot = nboot)
      out <- CorShrinkMatrix(cormat, nsamp=NULL, zscore_sd = zscore_sd,
                             thresh_up = thresh_up, thresh_down = thresh_down,
                             image = image,
                             tol=tol,
                             dosym=dosym,
                             image.control = image.control,
                             report_model = report_model,
                             maxiter = maxiter,
                             ash.control = ash.control)
    }
  }
  return(out)
}
