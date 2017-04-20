#' @title Adaptive shrinkage of a correlation matrix using ML estimation
#'
#' @description This function performs adaptive shrinkage of a sample correlation matrix using a
#' mixture normal prior on Fisher z-scores with fixed grid variances and the MLEs of
#' mixture proportions are calculated using the EM algorithm.
#'
#'
#' @param cormat The sample correlation matrix
#' @param nsamples The number of samples over which the correlation matrix is estimated.
#' @param sd_boot Boolean variable indicating if bootstrap sd to be used or not.
#' @param cor_transform_sd_vec If \code{sd_boot} is not NULL, the bootstrap standard errors
#' need to be provided. These standard errors can be obtained from \code{bootcorSE_calc()} function.
#' @param image if TRUE, plots an image of the shrunk and non-shrunk correlation matrices
#' @param tol The tolerance to check the difference between ash-cor only and ash-cor PD matrices.
#' @return Returns a shrunk version of the sample correlation matrix
#'         (the output is also a correlation matrix)
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216; doi: http://dx.doi.org/10.1101/038216
#' @keywords adaptive shrinkage, correlation
#' @importFrom reshape2 melt dcast
#' @import Matrix
#' @import ashr
#' @export


CorShrinkML <- function(cormat, nsamples, sd_boot = FALSE,
                        cor_transform_sd_vec = NULL,
                        thresh_up = 0.999, thresh_down = 0.001,
                        image=FALSE, tol=1e-06)
{
  cor_table <- reshape2::melt(cormat);
  cor_table_non_diag <- cor_table[which(cor_table[,1] != cor_table[,2]),];

  cor_table_non_diag.val <- cor_table_non_diag[,3];
  cor_table_non_diag.val[which(cor_table_non_diag.val >= thresh_up)]= thresh_up;
  cor_table_non_diag.val[which(cor_table_non_diag.val <= thresh_down)]= thresh_down;

  cor_transform_mean_vec=0.5*log((1+cor_table_non_diag.val)/(1-cor_table_non_diag.val))
  if(!sd_boot){
    cor_transform_sd_vec=rep(sqrt(1/(nsamples-3)), dim(cor_table_non_diag)[1]);
  }else{
    if(is.null(cor_transform_sd_vec)){
      stop("if sd_boot is not NULL, the user needs to provide a cor trasform sd vector")
    }
    cor_transform_sd_vec = cor_transform_sd_vec
  }
  options(warn=-1)
  fit=ashr::ash(cor_transform_mean_vec,cor_transform_sd_vec,
                mixcompdist="normal", optmethod="mixVBEM");

  ash_cor_vec=(exp(2*fit$result$PosteriorMean)-1)/(exp(2*fit$result$PosteriorMean)+1);

  newdata.table <- cor_table_non_diag;
  newdata.table[,3] <- ash_cor_vec;
  ash_cor_only <- reshape2::dcast(newdata.table, Var1~Var2, value.var = "value")[,-1];
  ash_cor_only[is.na(ash_cor_only)]=1;
  pd_completion <- Matrix::nearPD(as.matrix(ash_cor_only), conv.tol=tol);
  ash_cor_PD <- sweep(pd_completion$mat,diag(as.matrix(pd_completion$mat)), MARGIN=1,"/")
  if(image) {
    image(cormat)
    image(as.matrix(ash_cor_only))
  }
  if(all.equal(target=ash_cor_only, current=ash_cor_PD, tolerance=tol)==TRUE){
    cat("ash cor only and ash cor PD matrices are same")
  }else{
    cat("ash cor only and ash cor PD matrices are different")
  }
  ll <- list("ash_cor_only"= ash_cor_only, "ash_cor_PD"=ash_cor_PD)
  return(ll)
}

