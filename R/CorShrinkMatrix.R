#' @title Adaptive shrinkage of a matrix of pairwise correlations.
#'
#' @description This function performs adaptive shrinkage of a matrix of pairwise correlations
#' using a mixture normal prior on Fisher z-scores, with each component centered at the
#' same base level z-score value (0 for 0 base correlation) but a wide range of
#' data-driven component variances. The method is similar to the adaptive shrinkage method for
#' modeling false discovery rates proposed in Stephens 2016 (see reference).
#'
#'
#' @param cormat A matrix of pairwise correlations - not necessarily a correlation matrix.
#'               NAs in this matrix are treated as 0.
#' @param nsamp An integer or a matrix denoting the number of samples for
#'              each pair of variables over which the correlation has been computed.
#'              Only used when \code{zscore_sd} is not provided.
#' @param zscore_sd A matrix of the sandard error of the Fisher z-scores for each pair of
#'                 variables. May contain NA-s as well. The NA-s in this matrix must
#'                 match with the NAs in the \code{cormat} matrix. If provided, it is
#'                 used as default over the the asymptotic formulation using \code{nsamp}.
#'                 When set to NULL, asymptotic distribution of the Fisher z-scores is used
#'                 using \code{nsamp}.
#' @param thresh_up Upper threshold for correlations in \code{cormat}. Defaults to 0.99
#' @param thresh_down Lower threshold for correlations in \code{cormat}. Defaults to -0.99
#' @param image character. options for plotting the original or the corshrink matrix.
#'              If \code{image = "both"}, then the function outputs both the plot
#'              for original and shrunk correlationmatrix. If \code{image = "original"},
#'              then the function outputs the correlation plot for the original matrix only.
#'              If \code{image = "corshrink"}, then the function outputs the correlation plot
#'              for the CorShrink matrix only.If \code{image = "output"}, then the function
#'              outputs the saved ggplot figure without displaying it. If \code{image = "null"},
#'              no image is output. Defaults to "both".
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
#' @return If \code{report_model = FALSE}, returns a list with adaptively shrunk version
#'         of the sample correlation matrix both before (\code{ash_cor_only}) and after
#'         PD completion (\code{ash_cor_PD}). If \code{report_model = TRUE}, then the
#'         function also returns all the details of the adaptive shrinkage model output.
#'
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216;
#'              doi: http://dx.doi.org/10.1101/038216
#'
#' @examples
#'
#' data("pairwise_corr_matrix")
#' data("common_samples")
#' out <- CorShrinkMatrix(pairwise_corr_matrix, common_samples)
#'
#' @keywords adaptive shrinkage, correlation
#' @importFrom reshape2 melt dcast
#' @importFrom stats cor sd
#' @importFrom corrplot corrplot
#' @importFrom gridExtra grid.arrange
#' @importFrom utils modifyList
#' @importFrom Matrix nearPD
#' @importFrom ashr ash
#' @export


CorShrinkMatrix <- function(cormat, nsamp = NULL,
                        zscore_sd = NULL,
                        thresh_up = 0.99, thresh_down = - 0.99,
                        image = c("both", "original", "corshrink", "output", "null"),
                        tol=1e-06,
                        image.control = list(),
                        report_model = FALSE,
                        ash.control = list())
{

  if(length(image) > 1){
    image <- "null"
  }

  image.control.default <- list(method = "color", type =  "upper", add = FALSE,
                                col = colorRampPalette(c("blue", "white", "red"))(200), bg = "white",
                                title = "", is.corr = TRUE, diag = FALSE,
                                outline = FALSE, mar = c(0, 0, 0, 0), addgrid.col = NULL,
                                addCoef.col = NULL, addCoefasPercent = FALSE,
                                order = c("original", "AOE", "FPC", "hclust", "alphabet"),
                                hclust.method = c("complete", "ward",
                                                  "ward.D", "ward.D2", "single", "average",
                                                  "mcquitty", "median", "centroid"),
                                addrect = NULL, rect.col = "white", rect.lwd = 2, tl.pos = "td",
                                tl.cex = 0.8, tl.col = "black", tl.offset = 0.4, tl.srt = 90,
                                cl.pos = NULL, cl.lim = NULL, cl.length = NULL, cl.cex = 0.8,
                                cl.ratio = 0.15, cl.align.text = "c", cl.offset = 0.5, number.cex = 1,
                                number.font = 2, number.digits = NULL, addshade = c("negative",
                                "positive", "all"), shade.lwd = 1, shade.col = "white", p.mat = NULL,
                                sig.level = 0.05, insig = c("pch", "p-value", "blank", "n", "label_sig"),
                                pch = 4, pch.col = "black", pch.cex = 3, plotCI = c("n", "square",
                               "circle", "rect"), lowCI.mat = NULL, uppCI.mat = NULL, na.label = "?",
                                na.label.col = "white", win.asp = 1)

  image.control <- utils::modifyList(image.control.default, image.control)

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

  if(image == "original"){
     image_original <- TRUE
     image_corshrink <- FALSE
   }
  if(image == "corshrink"){
     image_corshrink <- TRUE
     image_original <- FALSE
  }
  if(image == "output"){
    image_original = TRUE
    image_corshrink = TRUE
  }
  if(image == "both"){
    par(mfrow=c(1,2))
    image_original = TRUE
    image_corshrink = TRUE
  }
  if(image == "null"){
    image_original = FALSE
    image_corshrink = FALSE
  }

   if(image_original) {
      do.call(corrplot::corrplot, append(list(corr = as.matrix(cormat)),
                                            image.control))
   }

    if(image_corshrink){
      do.call(corrplot::corrplot, append(list(corr = as.matrix(ash_cor_PD)),
                                             image.control))
    }

  ash_cor_PD <- as.matrix(ash_cor_PD)

  if(!is.null(rownames(cormat)) && !is.null(colnames(cormat))){
    rownames(ash_cor_only) <- rownames(cormat)
    colnames(ash_cor_only) <- colnames(cormat)
    rownames(ash_cor_PD) <- rownames(cormat)
    colnames(ash_cor_PD) <- colnames(cormat)
  }

  if(report_model){
    ll <- list("cor_before_PD"= ash_cor_only,
               "cor"= ash_cor_PD,
               "model" = fit)
  }else{
    ll <- list("cor_before_PD"= ash_cor_only, "cor"= ash_cor_PD)
  }

   # if(image == "both"){
   #   gridExtra::grid.arrange(out1, out2, nrow = 1)
   # }else if (image == "original"){
   #   print(out1)
   # }else if (image == "corshrink"){
   #   print(out2)
   # }else if (image == "output"){
   #   ll[["image"]] = out2
   # }

  return(ll)
}

