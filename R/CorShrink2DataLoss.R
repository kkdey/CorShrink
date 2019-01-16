#' @title CorShrink Loss method for correlation shrinkage
#' @description Performs the CorShrink-Loss approach of correlation shrinkage, that
#'              uses the Taylor series expansion of the empirical correlation matrix
#'              around the population correlation matrix to derive the distribution.
#'
#' @param data_with_missing The samples by features data matrix. May contain NA values.
#' @param alpha The tuning parameter for the L-1 scaling.
#'
#' @examples
#' data("sample_by_feature_data")
#' out = CorShrink2DataLoss(sample_by_feature_data, alpha = 1)
#' corrplot::corrplot(as.matrix(out), diag = FALSE,
#'         col = colorRampPalette(c("blue", "white", "red"))(200),
#'         tl.pos = "td", tl.cex = 0.4, tl.col = "black",
#'         rect.col = "white",na.label.col = "white",
#'         method = "color", type = "upper")
#' out2 = CorShrink2DataLoss(sample_by_feature_data, alpha = 10)
#' corrplot::corrplot(out2, diag = FALSE,
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


CorShrink2DataLoss <- function(data_with_missing, alpha= 1){

  library(CVXR)
  ##################  Building matrix of common samples for pairwise comparisons  ####################

  binary_indicator = matrix(1, nrow(data_with_missing), ncol(data_with_missing))
  binary_indicator[is.na(data_with_missing)]= 0
  common_samples = t(binary_indicator)%*%binary_indicator
  diag(common_samples) = 0

  #################  Pairwise correlations computation  ###############################

  pairwise_cor = cor(data_with_missing, use = "pairwise.complete.obs")
  pairwise_cor[is.na(pairwise_cor)] = 0
  pairwise_cor[pairwise_cor > 0.95] = 0.95
  pairwise_cor[pairwise_cor < -0.95] = -0.95
  diag(pairwise_cor) = 1
  common_samples[common_samples <= 2] = 2

  ################# Computing sample Fisher Z scores   ###########################

  pairwise_zscores = apply(pairwise_cor, c(1,2), function(x) return (0.5*log((1+x)/(1-x))))
  diag(pairwise_zscores) = 0

  ################  Samples adjustment - bias and sd   ########################

  samples_adjustment = 1/(common_samples - 1)  + 2/(common_samples - 1)^2
  bias_cor = samples_adjustment*pairwise_cor*(1 + pairwise_cor)*(1 - pairwise_cor)
  sd_cor = sqrt(samples_adjustment)*(1+pairwise_cor)*(1-pairwise_cor)
  estimate_cor = pairwise_cor + bias_cor
  inverse_sd_cor = 1/sd_cor
  diag(inverse_sd_cor) = 0

  R <- Semidef(dim(pairwise_cor)[1])
  obj <- Minimize(alpha*p_norm(R, 1) + p_norm(mul_elemwise(inverse_sd_cor, square(R - estimate_cor)), 1))
  constraints = list(diag(R) == 1)
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  R_hat = cov2cor(as.matrix(result$getValue(R)))
  return(R_hat)
}
