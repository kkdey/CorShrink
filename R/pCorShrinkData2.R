#' @title Adaptive shrinkage of a partial correlation matrix
#'
#' @description This function performs adaptive shrinkage of a matrix of
#' partial correlations using a similar model as \code{CorShrinkMatrix}, but
#' with the option of being aided by bagging on the data matrix.
#'
#' @param data The samples by features data matrix. May contain NA values.
#' @param type The type of bagging method used, if at all, for estimating the
#'             partial corelation matrix. The three possible types are -
#'             \code{no_bag}, \code{bag_cor} and \code{bag_pcor} - indicating
#'             if bagging is at all used, or if is performed at the level of
#'             correlations before taking partial correlation matroix on the
#'             bagged correlation values, or at the level of partial correlations
#'             directly.
#' @param nbags If \code{type} not equal to \code{no_bag}, how many bootstrap
#'              samples are generated to get the bagging estimate.
#' @param nboot The number of bootstrap samples used to compute the bootstrap
#'              standard errors of the partial correlations used for shrinking
#'              in \code{CorShrink}.
#' @param thresh_up Upper threshold for correlations in \code{cormat}. Defaults to 0.99
#' @param thresh_down Lower threshold for correlations in \code{cormat}. Defaults to -0.99
#' @param tol The tolerance chosen to check how far apart the CorShrink matrix is from the nearest
#'            positive definite matrix before applying PD completion.
#' @param report_model  if TRUE, outputs the full adaptive shrinkage output, else outputs the shrunken vector.
#'                      Defaults to FALSE.
#' @param ash.control The control parameters for adaptive shrinkage
#'
#' @return If \code{report_model = FALSE}, returns the adaptively shrunk
#'         partial correlation matrix. If \code{report_model = TRUE}, then the
#'         function also returns all the details of the adaptive shrinkage model
#'         output.
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216;
#'              doi: http://dx.doi.org/10.1101/038216
#'
#' @keywords adaptive shrinkage, correlation
#' @importFrom reshape2 melt dcast
#' @importFrom stats cor sd
#' @importFrom corrplot corrplot
#' @importFrom gridExtra grid.arrange
#' @importFrom utils modifyList
#' @importFrom Matrix nearPD
#' @importFrom ashr ash
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics par
#' @export


pCorShrinkData2 <- function(data,
                           type = c("no_bag", "bag_cor", "bag_pcor"),
                           nbags = 10, nboot = 50,
                           thresh_up = 0.99, thresh_down = - 0.99,
                           tol = 1e-06,
                           report_model = FALSE,
                           ash.control = list()){

  if(length(type) > 1 ){
    type <- "no_bag"
  }

  if(type == "bag_cor"){
    message("using bagging at correlation level")
  }else if (type == "bag_pcor"){
    message("using bagging at the level of partial
            correlations")
  }else{
    type = "no_bag"
    message("using no bagging")
  }

  if(type == "no_bag"){
    est_pcor <- cor2pcor1(data)
  }

  if(type == "bag_cor"){
    est_pcor <- cor2pcor3(data, num_bags = nbags)
  }

  if(type == "bag_pcor"){
    est_pcor <- cor2pcor2(data, num_bags = nbags)
  }

  pcor_bag_summ <- array(0, c(dim(est_pcor)[1], dim(est_pcor)[2], nboot))
  for(num in 1:nboot){
    data_samp <- data[sample(1:dim(data)[1], dim(data)[1], replace = TRUE), ]
    if(type == "no_bag"){
      pcor_bag <- cor2pcor1(data_samp)
    }else if(type == "bag_cor"){
      pcor_bag <- cor2pcor3(data_samp)
    }else if (type == "bag_pcor"){
      pcor_bag <- cor2pcor2(data_samp)
    }
    pcor_bag_summ[,,num] <- pcor_bag
  }

  tmp <- atanh(pcor_bag_summ)
  tmp[tmp==Inf] = NA
  sd_pcor_bag <- apply(tmp, c(1,2), function(x) return(sd(x)))

  out <- CorShrinkMatrix(est_pcor, zscore_sd = sd_pcor_bag,
                         tol = tol,
                         report_model = report_model,
                         ash.control = ash.control)

  return(out)
}

cor2pcor1 <- function(data){
  return(cor2pcor(cor(data)))
}

cor2pcor2 <- function(data, num_bags = 50){
  pS_bag <- 0
  for(num in 1:num_bags){
    data_samp <- data[sample(1:dim(data)[1], dim(data)[1], replace = TRUE), ]
    S_bag <- cor(data_samp)
    pS_bag <- pS_bag + corpcor::cor2pcor(S_bag)
  }
  pS_bag <- pS_bag/num_bags
  return(pS_bag)
}

cor2pcor3 <- function(data, num_bags = 50){
  S_bag <- 0
  for(num in 1:num_bags){
    data_samp <- data[sample(1:dim(data)[1], dim(data)[1], replace = TRUE), ]
    S_bag <- S_bag + cor(data_samp)
  }
  pS_bag <- cor2pcor((S_bag/num_bags))
  return(pS_bag)
}


