#' @title Generate a covariance matrix
#'
#' @description Wrapper for genPositiveDefMat function in clusterGeneration package to generate covariance matrix
#'
#' @param dim The dimensionality of the covariance matrix to generate
#' @return Returns a covariance matrix of a given dimensionality along with its eigenvalues vector.
#'
#' @keywords covariance matrix, eigenvalues
#' @export

generate_cov <- function(dim){
  pop_cov_class <- clusterGeneration::genPositiveDefMat(dim, covMethod="eigen")
  pop_cov_eigen <- pop_cov_class$egvalues;
  pop_cov_Sigma <- pop_cov_class$Sigma;
  ll <- list("Sigma"=pop_cov_Sigma, "eigen"=pop_cov_eigen);
  return(ll)
}

