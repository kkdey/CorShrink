#' @title Adaptive shrinkage of partial correlations from a data matrix
#'
#' @description Performs adaptive shrinkage of the sample inverse covariances and
#' consequently partial correlation matrix starting from a data matrix ( no NA
#' permitted unlike \code{CorShrinkData}). The procedure of shrinkage combines
#' inverse covariance derivations from the ISEE algorithm of [Fan and Lv, 2016]
#' with the \code{CorShrink} formulation.
#'
#' @param dat The samples by features data matrix. Does not permit NA entries in X.
#' @param permu_type Determines if the columns of the data matrix (the features)
#'                   are permuted randomly before blocking in the ISEE algorithm
#'                   [see reference] or not.Takes one of two character input values
#'                   - "random" and "ordered" - depending on if the columns are permuted
#'                   or not.
#' @param thresh_up Upper threshold for correlations. Defaults to 0.99
#' @param thresh_down Lower threshold for correlations. Defaults to -0.99.
#' @param ash.control The control parameters for adaptive shrinkage
#'
#' @return Returns an adaptively shrunk version of the inverse covariance matrix and
#'         the partial correlation matrix estimate.
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216;
#'              doi: http://dx.doi.org/10.1101/038216
#'              Fan, Y. and Lv, J., 2016. Innovated scalable efficient estimation
#'              in ultra-large Gaussian graphical models. The Annals of Statistics, 44(5), pp.2098-2126.
#' @keywords adaptive shrinkage, inverse covariance, partial correlation
#' @importFrom stats cor sd lm
#' @importFrom corpcor cor2pcor
#' @export
#'


pCorShrinkData <- function(dat,
                           permu_type = c("random", "ordered"),
                           thresh_up = 0.99,
                           thresh_down = 0.01,
                           ash.control = list()){

  if(length(permu_type) > 1){
    permu_type <- "ordered"
  }

  if(permu_type == "ordered"){
    permu <- 1:dim(dat)[2]
  }else{
    permu <- sample.int(dim(dat)[2])
  }

  dat <- dat[, permu]
  K <- 2
  dat_scaled = scale(dat, center=F, scale=T)/sqrt(nrow(dat)-1)
  y.norm = attr(dat_scaled,"scale")*sqrt(nrow(dat)-1)
  i.index <- seq(1,ncol(dat)-1,by=K)
  dat_tilde = matrix(0, nrow=nrow(dat), ncol=ncol(dat))
  p = ncol(dat)

  for (k in 1:length(i.index)){
    i = i.index[k]

    if(i==p-2){
      j=c(p-1,p)
    }else{
      j=i+1
    }

    block.ind = c(i,j)
    temp= lm(dat_scaled[, block.ind] ~ dat_scaled[, -block.ind] - 1)
    beta.coef = temp$coefficients

    temp2 = scale(dat_scaled[,-block.ind] %*% beta.coef, center=F, scale=1/y.norm[block.ind])
    epsi <- dat[, block.ind] - temp2
    omega <- solve(t(epsi)%*%epsi/n)
    dat_tilde[,block.ind] =  as.matrix(epsi%*%omega)
  }

  dat_tilde = dat_tilde[,order(permu)]
  out <- CorShrinkData(dat_tilde, image = "null", sd_boot = FALSE,
                       thresh_up = thresh_up,
                       thresh_down = thresh_down,
                       ash.control = ash.control)
  pcorShrink <- -as.matrix(out$cor)
  diag(pcorShrink) <- 1
  return(pcorShrink)
}

