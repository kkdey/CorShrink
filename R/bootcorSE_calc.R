#' @title Calculate bootstrap SE for the Fisher z-scores of correlations
#'
#' @description computes the bootstrap SE of the Fisher z-scores to be used in \code{CorShrink-ML}
#' or \code{CorShrink-VEM}.
#'
#' @param data The samples by features data matrix.
#' @param nboot The number of bootstrap samples.
#' @param thresh_up The upward threshold for correlations
#' @param thresh_down The downward threshold for correlations.
#'
#' @return Returns standard errors for fisher z-scores for the correlations.
#' @keywords adaptive shrinkage, correlation
#' @importFrom reshape2 melt dcast
#' @import Matrix
#' @import ashr
#' @export



bootcorSE_calc <- function(data, nboot = 50, thresh_up = 0.999, thresh_down = 0.001){

  boot_mat <- matrix(0, nboot, dim(data)[2]*(dim(data)[2]-1))

  for(num in 1:nboot){
    boot_indices <- sample(1:dim(data)[1], dim(data)[1], replace = TRUE)
    data_boot <- data[boot_indices,]
    cov_boot <- cov(data_boot)
    cor_sample_boot <- cov2cor(cov_boot)

    cor_table_boot <- reshape2::melt(cor_sample_boot);
    cor_table_non_diag_boot <- cor_table_boot[which(cor_table_boot[,1] != cor_table_boot[,2]),];

    cor_table_non_diag_boot.val <- cor_table_non_diag_boot[,3];
    cor_table_non_diag_boot.val[which(cor_table_non_diag_boot.val>= thresh_up)] = thresh_up
    cor_table_non_diag_boot.val[which(cor_table_non_diag_boot.val<= thresh_down)] = thresh_down

    cor_transform_mean_vec_boot=0.5*log((1+cor_table_non_diag_boot.val)/(1-cor_table_non_diag_boot.val))

    boot_mat[num, ] <- cor_transform_mean_vec_boot
    cat("Finished Bootstrap :", num, "\n")
  }

  sd_boot <- apply(boot_mat, 2, function(x) return(sd(x)))

  return(sd_boot)
}
