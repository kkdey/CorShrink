#' @title Box-Shrinkage of correlations from a data matrix with missing data
#' @description Performs a box shrinkage of correlations from a data matrix with missing data
#' using the Fisher Z-score formulation
#'
#' @param data_with_missing The samples by features data matrix. May contain NA values.
#'
#' @examples
#' data("sample_by_feature_data")
#' out = CorShrink2Data(sample_by_feature_data)
#' corrplot::corrplot(as.matrix(out), diag = FALSE,
#'         col = colorRampPalette(c("blue", "white", "red"))(200),
#'         tl.pos = "td", tl.cex = 0.4, tl.col = "black",
#'         rect.col = "white",na.label.col = "white",
#'         method = "color", type = "upper")
#'
#' @keywords box shrinkage, correlation
#' @import CVXR
#' @importFrom corrplot corrplot
#' @importFrom stats cor sd cov2cor
#' @export

CorShrink2Data <- function(data_with_missing){

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

  ###############  Bounds on the correlations  ####################################

  bound1 = 12*exp(2*pairwise_zscores)/((exp(2*pairwise_zscores) + 1)^2)
  zscores_sd = sqrt(1/(common_samples - 1) + 2/(common_samples - 1)^2)
  bound2 = bound1*zscores_sd
  bound3 = zscores_sd^2*2*sqrt(3)
  overall_bound = bound2 + bound3
  constrained_overall_bound = apply(overall_bound, c(1,2), function(x) return(min(2,x)))
  diag(constrained_overall_bound) = 0

  ###############  Convec optimization  ######################
  library(CVXR)
  R <- Semidef(dim(pairwise_cor)[1])
  obj <- Minimize(sum(abs(R)))
  constraints = list(diag(R) == 1, R <= pairwise_cor + constrained_overall_bound,
                     R >= pairwise_cor - constrained_overall_bound)
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  R_hat = as.matrix(result$getValue(R))
  R_hat = cov2cor(R_hat)
  return(R_hat)
}
