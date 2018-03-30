#' @title Adaptive shrinkage of a vector of correlations
#'
#' @description This function performs adaptive shrinkage of a
#' vector of sample correlations using a mixture normal prior on
#' Fisher z-scores , with each component centered at the
#' same base level z-score value (0 for 0 base correlation) but a
#' wide range of data-driven component variances. The method is similar to the
#' adaptive shrinkage method for modeling false discovery rates proposed in
#' Stephens 2016 (see reference).
#'
#' @param corvec A vector of sample correlations (may contain NAs)
#' @param nsamp_vec A vector of the number of samples over which correlations
#' for each cell of the vector are estimated.
#' @param zscore_sd_vec A vector of the sandard error of the Fisher z-scores for each cell
#'                 in the vector. May contain NA-s as well. The NA-s in this matrix must
#'                 match with the NAs in \code{corvec}. If provided, it is
#'                 used as default. When set to NULL, asymptotic distribution of the
#'                 Fisher z-scores is used using \code{nsamp_vec}.
#' @param thresh_up Upper threshold for correlations in \code{corvec}. Defaults to 0.99
#' @param thresh_down Lower threshold for correlations in \code{corvec}. Defaults to -0.99
#'
#' @param report_model  if TRUE, outputs the full adaptive shrinkage output, else outputs
#'                      the shrunken vector. Defaults to FALSE.
#' @param ash.control The control parameters for adaptive shrinkage
#'
#' @return If \code{report_model = FALSE}, returns an adaptively shrunk version
#'         of the vector of correlations. If \code{report_model = TRUE}, then the
#'         function also returns all the details of the adaptive shrinkage model output.
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216;
#'              doi: http://dx.doi.org/10.1101/038216
#' @examples
#'
#'  cor_vec <- c(-0.56, -0.4, 0.02, 0.2, 0.9, 0.8, 0.3, 0.1, 0.4)
#'  nsamp_vec <- c(10, 20, 30, 4, 50, 60, 20, 10, 3)
#'  out <- CorShrinkVector(corvec = cor_vec, nsamp_vec = nsamp_vec)
#'
#' @keywords adaptive shrinkage, correlation
#' @importFrom ashr ash
#' @importFrom stats cor sd
#' @importFrom utils modifyList
#' @export
CorShrinkVector <- function (corvec, nsamp_vec,
                             zscore_sd_vec = NULL,
                             thresh_up = 0.99, thresh_down = - 0.99,
                             report_model = FALSE,
                             ash.control = list()){


  if(!is.numeric(corvec)){
    stop("corvec must be a numeric vector")
  }

  if(is.null(nsamp_vec) && is.null(zscore_sd_vec)){
    stop("One of the two vectors: nsamp_vec  or zscore_sd_vec must be non NULL")
  }

  if(!is.null(nsamp_vec) && !is.null(zscore_sd_vec)){
    warning("Both nsamp_vec and zscore_sd_vec are non NULL : switching to nsamp_vec")
    zscore_sd_vec <- NULL
  }

  if(!is.null(nsamp_vec) && !is.numeric(nsamp_vec)){
    stop("nsamp_vec is specified but is not a numeric vector")
  }

  if(!is.null(zscore_sd_vec) && !is.numeric(zscore_sd_vec)){
    stop("zscore_sd_vec is specified but is not a numeric vector")
  }

  if(max(corvec) > 1 | min(corvec) < -1){
    stop("some values in corvec lie outside the [-1, 1] range for correlations")
  }

  ash.control.default = list(pointmass = TRUE,
                             mixcompdist = "normal", nullweight = 10,
                             fixg = FALSE, mode = 0, optmethod = "mixEM",
                             prior = "nullbiased", gridmult = sqrt(2),
                             outputlevel = 2, alpha = 0,
                             df = NULL, control=list(K = 1, method=3,
                             square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1,
                             objfn.inc=1,tol=1.e-05, maxiter=100, trace=FALSE))

  ash.control <- utils::modifyList(ash.control.default, ash.control)

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
                               sebetahat = corvec_trans_sd),
                               ash.control))
 ash_corvec=(exp(2*fit$result$PosteriorMean)-1)/(exp(2*fit$result$PosteriorMean)+1);

 if(report_model){
   ll <- list("estimate" = ash_corvec,
              "model" = fit)
   return(ll)
 }else{
   return(ash_corvec)
 }
}
