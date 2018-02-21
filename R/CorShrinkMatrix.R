#' @title Adaptive shrinkage of a matrix of correlations.
#'
#' @description This function performs adaptive shrinkage of a sample correlation matrix using a
#' mixture normal prior on Fisher z-scores with wide range of grid variances.
#' The method is similar to the adaptive shrinkage method for modeling false discovery rates proposed
#' in Stephens 2016 (see reference).
#'
#'
#' @param cormat A table of correlations - not necessarily a correlation matrix.
#'               May contain NAs as well.
#' @param nsamp An integer or a matrix denoting the number of samples for
#'              each pair of variables over which the correlation has been computed.
#'              If the user specifies \code{zscore_sd}, then \code{nsamp} is set
#'              to NULL and is no longer used.
#' @param zscore_sd A matrix of the sandard error of the Fisher z-scores for each pair of
#'                 variables. May contain NA-s as well. The NA-s in this matrix must
#'                 match with the NAs in the \code{cormat} matrix. If provided, it is
#'                 used as default over the the asymptotic formulation using \code{nsamp}.
#'                 When set to NULL, asymptotic distribution of the Fisher z-scores is used
#'                 using \code{nsamp}.
#' @param thresh_up Upper threshold for correlations. Defaults to 0.99
#' @param thresh_down Lower threshold for correlations. Defaults to -0.99
#' @param image_original if TRUE, plots an image of the non-shrunk original matrix
#'                       of correlations.
#' @param image_corshrink if TRUE, plots an image of the shrunk matrix
#'                       of correlations.
#' @param tol The tolerance chosen to check how far apart the CorShrink matrix is from the nearest
#'            positive definite matrix before applying PD completion.
#'
#' @param image.control Control parameters for the image when
#'                      \code{image_original = TRUE} and/or \code{image_corshrink = TRUE}.
#'
#' @param report_model  if TRUE, outputs the full adaptive shrinkage output, else outputs the shrunken vector.
#'                      Defaults to FALSE.
#' @param ash.control The control parameters for adaptive shrinkage
#'
#' @return Returns a adaptively shrunk version of the sample correlation matrix before
#'         and after PD completion.
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216; doi: http://dx.doi.org/10.1101/038216
#'
#' @examples
#'
#' data("pairwise_corr_matrix")
#' data("common_samples")
#' out <- CorShrinkMatrix(pairwise_corr_matrix, common_samples, image_corshrink  = TRUE)
#'
#' @keywords adaptive shrinkage, correlation
#' @importFrom reshape2 melt dcast
#' @importFrom grDevices rgb
#' @importFrom graphics axis
#' @importFrom stats cor sd
#' @importFrom utils modifyList
#' @import Matrix
#' @import ashr
#' @import SQUAREM
#' @export


CorShrinkMatrix <- function(cormat, nsamp = NULL,
                        zscore_sd = NULL,
                        thresh_up = 0.99, thresh_down = - 0.99,
                        image_original=FALSE, image_corshrink = FALSE,
                        tol=1e-06,
                        image.control = list(),
                        report_model = FALSE,
                        ash.control = list())
{
  image.control.default <- list(x.las = 2,
                                x.cex = 0.7,
                                y.las = 2,
                                y.cex = 0.7,
                                main_original = "sample corr matrix",
                                main_corshrink = "CorShrink matrix",
                                cex.main = 1,
                                col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
                                      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000))))
  image.control <- modifyList(image.control.default, image.control)

  if(is.null(zscore_sd) && is.null(nsamp)){
    stop("User must provide wither an nsamp or a zscore_sd arguments")
  }

  cormat[is.na(cormat)] = 0

  if(is.null(rownames(cormat))){
    rownames(cormat) <- 1:dim(cormat)[1]
  }

  if(is.null(colnames(cormat))){
    colnames(cormat) <- 1:dim(cormat)[2]
  }

  if(!is.null(zscore_sd)){
    zscore_sd[is.na(zscore_sd)] = 0
  }

  if(!is.null(nsamp)){
    nsamp[is.na(nsamp)] <- 0
  }

  ##############  Set control parameters for adaptive shrinkage  ################

  ash.control.default = list(pointmass = TRUE,
                             mixcompdist = "normal", nullweight = 10,
                             fixg = FALSE, mode = 0, optmethod = "mixEM",
                             prior = "nullbiased", gridmult = sqrt(2),
                             outputlevel = 2, alpha = 0,
                             df = NULL, control = list(K = 1,
                             method=3, square=TRUE, step.min0=1, step.max0=1,
                             mstep=4, kr=1, objfn.inc=1,tol=1.e-05, maxiter=100, trace=FALSE))
  ash.control <- modifyList(ash.control.default, ash.control)

  ##################   vectorise the correlation matrix  ###################

  cor_table <- reshape2::melt(cormat);
  cor_table_non_diag <- cor_table[which(cor_table[,1] != cor_table[,2]),];

  ##################  thresholding very low or very high correlations ##########

  cor_table_non_diag.val <- cor_table_non_diag[,3];
  cor_table_non_diag.val[which(cor_table_non_diag.val >= thresh_up)]= thresh_up;
  cor_table_non_diag.val[which(cor_table_non_diag.val <= thresh_down)]= thresh_down;

  ##################  Compute Fisher z-transform  ###########################

  cor_transform_mean_vec=0.5*log((1+cor_table_non_diag.val)/(1-cor_table_non_diag.val))

  ################ Compute standard errors of Fisher z-transform  ############


  if(is.null(zscore_sd) && (length(nsamp) == 1)){
    nsamples <- as.numeric(round(nsamp))
    if(nsamples <= 2){
      stop("the number of samples <=2 for all cells, will result in Identity correlation matrix in CorShrink")
    }
    cor_transform_sd_vec=rep(sqrt(1/(nsamples-1) + 2/(nsamples - 1)^2), dim(cor_table_non_diag)[1]);
  }else if(is.null(zscore_sd) && is.matrix(nsamp)){

    nsamp_tab <- reshape2::melt(nsamp)
    nsamp_tab_non_diag <- nsamp_tab[which(nsamp_tab[,1] != nsamp_tab[,2]),];
    nsamp_vec <- nsamp_tab_non_diag[,3]
    index_zeros <- which(nsamp_vec <= 2)
    cor_transform_mean_vec[index_zeros] = 0;
    nsamp_vec_2 <- nsamp_vec
    nsamp_vec_2[index_zeros] <- 1.00001
    cor_transform_sd_vec <- sqrt(1/(nsamp_vec_2-1) + 2/(nsamp_vec_2 - 1)^2);

  }else{
    if(is.null(zscore_sd)){
      stop("if sd_boot is not NULL, the user needs to provide a cor trasform sd vector")
    }
    cor_transform_sd_mat <- reshape2::melt(as.matrix(zscore_sd))
    cor_transform_sd_non_diag <- cor_transform_sd_mat[which(cor_transform_sd_mat[,1] != cor_transform_sd_mat[,2]),];
    cor_transform_sd_vec <- cor_transform_sd_non_diag[,3]
    index_zeros <- which(cor_transform_sd_vec  == 0)
    cor_transform_sd_vec[index_zeros] <- 10^8
  }

  options(warn=-1)

  ##################   Adaptive Shrinkage (Fisher Z scores) ################

  fit=do.call(ashr::ash, append(list(betahat = cor_transform_mean_vec,
                                     sebetahat = cor_transform_sd_vec),
                                ash.control))


  ###############   Inverse Fisher z-score transformation  ##################

  ash_cor_vec=(exp(2*fit$result$PosteriorMean)-1)/(exp(2*fit$result$PosteriorMean)+1);

  newdata.table <- cor_table_non_diag;
  newdata.table[,3] <- ash_cor_vec;
  ash_cor_only <- reshape2::dcast(newdata.table, Var1~Var2, value.var = "value")[,-1];
  ash_cor_only[is.na(ash_cor_only)]=1;
  rownames(ash_cor_only) <- rownames(cormat)
  colnames(ash_cor_only) <- colnames(cormat)


  ###############  Positive definite matrix completion of corShrink #############

  pd_completion <- Matrix::nearPD(as.matrix(ash_cor_only), conv.tol=tol);
  ash_cor_PD <- sweep(pd_completion$mat,diag(as.matrix(pd_completion$mat)), MARGIN=1,"/")


  if(is.null(rownames(cormat))){
    rownames(cormat) <- 1:dim(cormat)[1]
  }
  if(is.null(colnames(cormat))){
    colnames(cormat) <- 1:dim(cormat)[2]
  }

  row_labs <- rownames(cormat)
  col_labs <- colnames(cormat)

   if(image_original) {
      image(as.matrix(cormat), col=image.control$col, main = image.control$main_original,
            cex.main=image.control$cex.main, xaxt = "n", yaxt = "n", zlim=c(-1,1))
      axis(1, at = seq(0, 1, length.out = ncol(cormat)),
           labels = row_labs, las=image.control$x.las, cex.axis = image.control$x.cex)
      axis(2, at = seq(0, 1, length.out = ncol(cormat)),
           labels = col_labs, las=image.control$y.las, cex.axis = image.control$y.cex)
   }

    if(image_corshrink){
      image(as.matrix(ash_cor_PD), col=image.control$col, main=image.control$main_corshrink,
            cex.main=image.control$cex.main, xaxt = "n", yaxt = "n", zlim=c(-1,1))
      axis(1, at = seq(0, 1, length.out = ncol(cormat)),
           labels = row_labs, las=image.control$x.las, cex.axis = image.control$x.cex)
      axis(2, at = seq(0, 1, length.out = ncol(cormat)),
           labels = col_labs, las=image.control$y.las, cex.axis = image.control$y.cex)
    }

  if(all.equal(target=ash_cor_only, current=ash_cor_PD, tolerance=tol)==TRUE){
    message("ash cor only and ash cor PD matrices are same")
  }else{
    message("ash cor only and ash cor PD matrices are different")
  }

  if(report_model){
    ll <- list("ash_cor_only"= ash_cor_only,
               "ash_cor_PD"=as.matrix(ash_cor_PD),
               "model" = fit)
    return(ll)
  }else{
    ll <- list("ash_cor_only"= ash_cor_only, "ash_cor_PD"=as.matrix(ash_cor_PD))
    return(ll)
  }
}

