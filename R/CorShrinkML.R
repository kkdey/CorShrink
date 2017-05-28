#' @title Adaptive shrinkage of a correlation matrix using ML estimation
#'
#' @description This function performs adaptive shrinkage of a sample correlation matrix using a
#' mixture normal prior on Fisher z-scores with fixed grid variances and the MLEs of
#' mixture proportions are calculated using the EM algorithm.
#'
#'
#' @param cormat The sample correlation matrix
#' @param nsamp_mat The number of samples/ matrix of samples over which the correlation matrix is estimated.
#' @param sd_boot Boolean variable indicating if bootstrap sd to be used or not.
#' @param cor_transform_sd_vec If \code{sd_boot} is not NULL, the bootstrap standard errors
#' need to be provided. These standard errors can be obtained from \code{bootcorSE_calc()} function.
#' @param image if TRUE, plots an image of the shrunk and non-shrunk correlation matrices
#' @param tol The tolerance to check the difference between ash-cor only and ash-cor PD matrices.
#' @param ash.control The control parameters for adaptive shrinkage
#'
#' @return Returns a shrunk version of the sample correlation matrix
#'         (the output is also a correlation matrix)
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216; doi: http://dx.doi.org/10.1101/038216
#' @keywords adaptive shrinkage, correlation
#' @importFrom reshape2 melt dcast
#' @import Matrix
#' @import ashr
#' @export


CorShrinkML <- function(cormat, nsamp_mat, sd_boot = FALSE,
                        cor_transform_sd_vec = NULL,
                        thresh_up = 0.99, thresh_down = - 0.99,
                        image=FALSE, tol=1e-06,
                        ash.control = list())
{
  cormat[is.na(cormat)] = 0
  ash.control.default = list(pointmass = TRUE, prior = "nullbiased", gridmult = 2,
                             mixcompdist = "normal", nullweight = 10,
                             outputlevel = 2, fixg = FALSE, optmethod="mixEM")
  ash.control <- modifyList(ash.control.default, ash.control)

  cor_table <- reshape2::melt(cormat);
  cor_table_non_diag <- cor_table[which(cor_table[,1] != cor_table[,2]),];

  cor_table_non_diag.val <- cor_table_non_diag[,3];
  cor_table_non_diag.val[which(cor_table_non_diag.val >= thresh_up)]= thresh_up;
  cor_table_non_diag.val[which(cor_table_non_diag.val <= thresh_down)]= thresh_down;

  cor_transform_mean_vec=0.5*log((1+cor_table_non_diag.val)/(1-cor_table_non_diag.val))

  if(!sd_boot && !is.matrix(nsamp_mat)){
    nsamples <- as.numeric(nsamp_mat)
    if(nsamples <= 2){
      stop("the number of samples <=2 for all cells, will result in 0 correlation matrix")
    }
    cor_transform_sd_vec=rep(sqrt(1/(nsamples-1) + 2/(nsamples - 1)^2), dim(cor_table_non_diag)[1]);
  }else if(!sd_boot && is.matrix(nsamp_mat)){

    nsamp_tab <- reshape2::melt(nsamp_mat)
    nsamp_tab_non_diag <- nsamp_tab[which(nsamp_tab[,1] != nsamp_tab[,2]),];
    nsamp_vec <- nsamp_tab_non_diag[,3]
    index_zeros <- which(nsamp_vec <= 2)
    cor_transform_mean_vec[index_zeros] = 0;
    nsamp_vec_2 <- nsamp_vec
    nsamp_vec_2[index_zeros] <- 1.0001
    cor_transform_sd_vec <- sqrt(1/(nsamp_vec_2-1) + 2/(nsamp_vec_2 - 1)^2);

  }else{
    if(is.null(cor_transform_sd_vec)){
      stop("if sd_boot is not NULL, the user needs to provide a cor trasform sd vector")
    }
    cor_transform_sd_vec = cor_transform_sd_vec
  }

  options(warn=-1)
  fit=do.call(ashr::ash, append(list(betahat = cor_transform_mean_vec,
                                     sebetahat = cor_transform_sd_vec),
                                ash.control))

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

