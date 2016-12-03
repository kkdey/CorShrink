#' @title Performs Shafer-Strimmer shrinkage on covariance matrix
#'
#' @description This function performs Shafer Strimmer shrinkage on a covariance matrix
#'
#' @param data The samples by features data matrix.
#' @return Returns a list with two items - the Shafer-Strimmer shrunk covariance matrix and the shrinkage intensity.
#'
#' @keywords Shafer-Strimmer shrinkage
#' @export

shafer_strimmer_shrinker <- function(data){
  nsamples <- dim(data)[1];
  nfeat <- dim(data)[2];
  colmean <- colMeans(data);
  covmat <- cov(data);
  cormat <- cov2cor(covmat);
  w <- array(0, c(nsamples, nfeat, nfeat));
  for(n in 1:nsamples){
    w[n,,] <- (data[n,] - colmean)%*% t(data[n,] - colmean);
  }
  var_hat_s <- (nsamples/(nsamples-1)^2) * apply(w, c(2,3), var);
  sum_var_hat_s <- sum(var_hat_s[row(var_hat_s)!=col(var_hat_s)])
  square_cor <- covmat^2;
  sum_s_square <- sum(square_cor[row(square_cor)!=col(square_cor)]);
  shrink_intensity <- sum_var_hat_s/sum_s_square;

  if(shrink_intensity > 1){
    shafer_shrink <- diag(diag(covmat));
  }
  if(shrink_intensity < 0){
    shafer_shrink <- covmat;
  }
  else{
    shafer_shrink <- diag(diag(covmat))
    shafer_shrink[row(shafer_shrink)!=col(shafer_shrink)] <- (1- shrink_intensity)*covmat[row(covmat)!=col(covmat)]
  }
  ll <- list("mat"=shafer_shrink, "shrink_intensity"=shrink_intensity)
  return(ll)
}
