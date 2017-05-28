#' @title Adaptive shrinkage of a vector of correlations using ML estimation
#'
#' @description This function performs adaptive shrinkage of a sample correlation matrix using a
#' mixture normal prior on Fisher z-scores with fixed grid variances and the MLEs of
#' mixture proportions are calculated using the EM algorithm.
#'
#'
#' @param corvec The vector of sample correlations
#' @param nsamp_vec vector of the number of samples over which the correlation matrix is estimated.
#' @param sd_boot Boolean variable indicating if bootstrap sd to be used or not.
#' @param cor_transform_sd_vec If \code{sd_boot} is not NULL, the vector of bootstrap standard errors
#' need to be provided. These standard errors can be obtained from \code{bootcorSE_calc()} function.
#' @param thresh_up The upper threshold for correlations
#' @param thresh_down The lower threshold for correlations
#' @param ash.control The control parameters for adaptive shrinkage
#'
#' @return Returns an adaptively shrunk version of the vector of correlations.
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216; doi: http://dx.doi.org/10.1101/038216
#' @keywords adaptive shrinkage, correlation
#' @import ashr
#' @export


CorShrinkMLvec <- function (corvec, nsamp_vec, sd_boot = FALSE,
                            cor_transform_sd_vec = NULL,
                            thresh_up = 0.99, thresh_down = - 0.99,
                            ash.control = list()){

  ash.control.default = list(pointmass = TRUE, prior = "nullbiased", gridmult = 2,
                             mixcompdist = "normal", nullweight = 10,
                             outputlevel = 2, fixg = FALSE, mode =0,
                             optmethod="mixEM")

  ash.control <- modifyList(ash.control.default, ash.control)

  corvec[corvec >= thresh_up] <- thresh_up
  corvec[corvec <= thresh_down] <- thresh_down

  corvec_trans <- 0.5* log ((1+corvec)/(1-corvec))

  if(!sd_boot && (length(nsamp_vec)==1)){
    nsamples <- as.numeric(nsamp_vec)
    if(nsamples <= 2){
      ash_corvec <- rep(0, length(corvec))
      return(ash_corvec)
    }
    corvec_trans_sd =rep(sqrt(1/(nsamples-1) + 2/(nsamples - 1)^2), length(corvec_trans));
  }else if(!sd_boot && is.vector(nsamp_vec)){
    index_zeros <- which(nsamp_vec <= 2)
    nsamp_vec_2 <- nsamp_vec
    nsamp_vec_2[index_zeros] <- 1.0001
    corvec_trans[index_zeros] <- 0
    corvec_trans_sd =sqrt(1/(nsamp_vec_2-1) + 2/(nsamp_vec_2 - 1)^2);
  }else{
    if(is.null(cor_transform_sd_vec)){
      stop("if sd_boot is not NULL, the user needs to provide a cor trasform sd vector")
    }
    corvec_trans_sd = cor_transform_sd_vec
 }

 options(warn=-1)
 fit=do.call(ashr::ash, append(list(betahat = corvec_trans,
                               sebetahat = corvec_trans_sd),
                               ash.control))
 ash_corvec=(exp(2*fit$result$PosteriorMean)-1)/(exp(2*fit$result$PosteriorMean)+1);

 return(ash_corvec)
}
