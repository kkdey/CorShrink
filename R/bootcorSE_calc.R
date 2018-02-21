#' @title Calculate bootstrap SE for the Fisher z-scores of correlations
#'
#' @description computes the bootstrap SE of the Fisher z-scores to be used in \code{CorShrink-ML}
#' or \code{CorShrink-VEM}.
#'
#' @param data The samples by features data matrix.
#' @param nboot The number of bootstrap samples.
#' @param thresh_up The upward threshold for correlations
#' @param thresh_down The downward threshold for correlations.
#' @param verbose To print the status of Bootstrap runs
#'
#' @return Returns standard errors for fisher z-scores for the correlations.
#' @keywords adaptive shrinkage, correlation
#' @examples
#'
#' data("sample_by_feature_data")
#' zscoreSDmat <- bootcorSE_calc(sample_by_feature_data, verbose = FALSE)
#'
#' @importFrom reshape2 melt dcast
#' @export



bootcorSE_calc <- function(data, nboot = 50, thresh_up = 0.999,
                           thresh_down = 0.001, verbose = TRUE){

  boot_mat <- matrix(0, nboot, dim(data)[2]*(dim(data)[2]-1))

  for(num in 1:nboot){
    boot_indices <- sample(1:dim(data)[1], dim(data)[1], replace = TRUE)
    data_boot <- data[boot_indices,]
    cor_sample_boot <- cor(data_boot, use = "pairwise.complete.obs")

    cor_table_boot <- reshape2::melt(cor_sample_boot);
    cor_table_non_diag_boot <- cor_table_boot[which(cor_table_boot[,1] != cor_table_boot[,2]),];

    cor_table_non_diag_boot.val <- cor_table_non_diag_boot[,3];
    cor_table_non_diag_boot.val[which(cor_table_non_diag_boot.val>= thresh_up)] = thresh_up
    cor_table_non_diag_boot.val[which(cor_table_non_diag_boot.val<= thresh_down)] = thresh_down

    cor_transform_mean_vec_boot=0.5*log((1+cor_table_non_diag_boot.val)/(1-cor_table_non_diag_boot.val))

    boot_mat[num, ] <- cor_transform_mean_vec_boot
    if(verbose){
      cat("Finished Bootstrap :", num, "\n")
    }
  }

  sd_boot <- apply(boot_mat, 2, function(x) return(sd(x, na.rm=TRUE)))
  dt <- cbind.data.frame(cor_table_non_diag_boot[,1:2], sd_boot)
  temp <- reshape2::dcast(dt, Var1 ~ Var2)
  sd_boot_mat <- temp[,-1]
  rownames(sd_boot_mat)  <- temp[,1]
  diag(sd_boot_mat) <- 1
  sd_boot_mat[is.na(sd_boot_mat)] <- 10^7
  return(as.matrix(sd_boot_mat))
}
