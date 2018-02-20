#' @title Adaptive shrinkage of a vector of correlations
#'
#' @description This function performs adaptive shrinkage of a
#' vector of sample correlations using a mixture normal prior on
#' Fisher z-scores with a wide range of grid variances.
#'
#'
#' @param corvec A vector of sample correlations (may contain NAs)
#' @param nsamp_vec A vector of the number of samples over which correlations
#' for each cell of the vector are estimated.
#' @param zscore_sd_vec A vector of the sandard error of the Fisher z-scores for each cell
#'                 in the vector. May contain NA-s as well. The NA-s in this matrix must
#'                 match with the NAs in the \code{corvec} matrix. If provided, it is
#'                 used as default. When set to NULL, asymptotic distribution of the Fisher z-scores is used
#'                 using \code{nsamp_vec}.
#' @param thresh_up Upper threshold for correlations. Defaults to 0.99
#' @param thresh_down Lower threshold for correlations. Defaults to -0.99
#' @param optmethod The optimization method for EM algorithm - can be one of
#'                  two techniques \code{mixEM} (mixture EM) and \code{mixVBEM}
#'                  (mixture Variational Bayes EM) approaches.The default approach
#'                  is \code{mixEM}.
#' @param ash.control The control parameters for adaptive shrinkage
#'
#' @return Returns an adaptively shrunk version of the vector of correlations.
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216; doi: http://dx.doi.org/10.1101/038216
#' @examples
#'
#'  cor_vec <- c(-0.56, -0.4, 0.02, 0.2, 0.9, 0.8, 0.3, 0.1, 0.4)
#'  nsamp_vec <- c(10, 20, 30, 4, 50, 60, 20, 10, 3)
#'  out <- CorShrinkVector(corvec = cor_vec, nsamp_vec = nsamp_vec, optmethod = "mixEM")
#'
#' @keywords adaptive shrinkage, correlation
#' @import ashr
#' @import SQUAREM
#' @importFrom stats cor sd
#' @importFrom utils modifyList
#' @export


CorShrinkVector <- function (corvec, nsamp_vec,
                             zscore_sd_vec = NULL,
                             thresh_up = 0.99, thresh_down = - 0.99,
                             optmethod = "mixEM",
                             ash.control = list()){

  ash.control.default = list(pointmass = TRUE,
                             mixcompdist = "normal", nullweight = 10,
                             fixg = FALSE, mode = 0,
                             prior = "nullbiased", gridmult = sqrt(2),
                             outputlevel = 2, alpha = 0,
                             df = NULL, control=list(K = 1, method=3,
                             square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1,
                             objfn.inc=1,tol=1.e-05, maxiter=100, trace=FALSE))

  ash.control <- modifyList(ash.control.default, ash.control)

  corvec[is.na(corvec)] = 0

  if(!is.null(zscore_sd_vec)){
    zscore_sd_vec[is.na(zscore_sd_vec)] = 0
  }

  if(!is.null(nsamp_vec)){
    nsamp_vec[is.na(nsamp_vec)] <- 0
  }


  corvec[corvec >= thresh_up] <- thresh_up
  corvec[corvec <= thresh_down] <- thresh_down

  corvec_trans <- 0.5* log ((1+corvec)/(1-corvec))

  if(is.null(zscore_sd_vec) && (length(nsamp_vec)==1)){
    nsamples <- as.numeric(nsamp_vec)
    if(nsamples <= 2){
      ash_corvec <- rep(0, length(corvec))
      return(ash_corvec)
    }
    corvec_trans_sd =rep(sqrt(1/(nsamples-1) + 2/(nsamples - 1)^2), length(corvec_trans));
  }else if(is.null(zscore_sd_vec) && is.vector(nsamp_vec)){
    index_zeros <- which(nsamp_vec <= 2)
    nsamp_vec_2 <- nsamp_vec
    nsamp_vec_2[index_zeros] <- 1.0001
    corvec_trans[index_zeros] <- 0
    corvec_trans_sd =sqrt(1/(nsamp_vec_2-1) + 2/(nsamp_vec_2 - 1)^2);
  }else{
    if(is.null(zscore_sd_vec)){
      stop("the user needs to provide a zscore sd vector")
    }
    corvec_trans_sd = zscore_sd_vec
 }

 options(warn=-1)
 fit=do.call(ashr::ash, append(list(betahat = corvec_trans,
                               sebetahat = corvec_trans_sd,
                               optmethod = optmethod),
                               ash.control))
 ash_corvec=(exp(2*fit$result$PosteriorMean)-1)/(exp(2*fit$result$PosteriorMean)+1);

 return(ash_corvec)
}
