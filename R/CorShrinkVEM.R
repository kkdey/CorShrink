#' @title Adaptive shrinkage of a correlation matrix with Variational EM on the mixture proportions
#'
#' @description This function performs adaptive shrinkage of a sample correlation matrix using mixture normal
#' prior with known grid variances and a Dirichlet distribution model for the mixture proportions.
#' with high weights on the null component. The model is estimated through Variational EM on the mixture porportions.
#'
#' @param cormat The sample correlation matrix
#' @param nsamples The number of samples over which the correlation matrix is estimated.
#' @param image if TRUE, plots an image of the shrunk and non-shrunk correlation matrices
#' @param tol The tolerance to check the difference between ash-cor only and ash-cor PD matrices.
#' @param nullweight The weight of the null component of the ash prior.
#' @param null.comp The number of null components. Default is 1.
#'
#' @return Returns a shrunk version of the sample correlation matrix (the output is also a correlation matrix)
#' @references  False Discovery Rates: A New Deal. Matthew Stephens bioRxiv 038216; doi: http://dx.doi.org/10.1101/038216
#' @keywords adaptive shrinkage, correlation
#' @importFrom reshape2 melt dcast
#' @import Matrix
#' @import ashr
#' @export

CorShrinkVEM <- function(cormat, nsamples, image=FALSE, tol=1e-06,
                     nullweight = 10, null.comp =1){

  cor_table <- reshape2::melt(cormat);
  cor_table_non_diag <- cor_table[which(cor_table[,1] !=cor_table[,2]),];
  cor_table_non_diag.val <- cor_table_non_diag[,3];
  cor_table_non_diag.val[which(cor_table_non_diag.val==1)]=(1- 1e-7);
  cor_transform_mean_vec=0.5*log((1+cor_table_non_diag.val)/(1-cor_table_non_diag.val))
  cor_transform_sd_vec=rep(sqrt(1/(nsamples-3)), dim(cor_table_non_diag)[1]);
  options(warn=-1)

  betahat <- cor_transform_mean_vec
  sebetahat <- cor_transform_sd_vec
  completeobs <- which(!is.na(betahat));
  gridmult = sqrt(2);
  mixsd = autoselect.mixsd(betahat[completeobs]-mean(betahat[completeobs]),sebetahat[completeobs],gridmult)
  mixsd = c(0,mixsd)
  null.comp = which.min(mixsd) #which component is the "null"

  k = length(mixsd)
  prior = "nullbiased";
  prior = setprior(prior,k,nullweight,null.comp)
  randomstart = FALSE
  n = length(betahat)
  pi = initpi(k,n,null.comp,randomstart)

  g=normalmix(pi,rep(0,k),mixsd)
  df <- NULL
  control=list()
  pi.fit=estimate_mixprop(betahat[completeobs],sebetahat[completeobs],
                          g,prior,null.comp=null.comp,
                          optmethod="mixVBEM",df=df,control=control)

  postMean.fit <- postmean(pi.fit$g,betahat[completeobs],sebetahat[completeobs], v=NULL)

  ash_cor_vec=(exp(2*postMean.fit)-1)/(exp(2*postMean.fit)+1);

  newdata.table <- cor_table_non_diag;
  newdata.table[,3] <- ash_cor_vec;
  ash_cor_only <- reshape2::dcast(newdata.table, Var1~Var2, value.var = "value")[,-1];
  ash_cor_only[is.na(ash_cor_only)]=1;
  pd_completion <- Matrix::nearPD(as.matrix(ash_cor_only), conv.tol=1e-06);
  ash_cor_PD2 <- sweep(pd_completion$mat,diag(as.matrix(pd_completion$mat)), MARGIN=1,"/")

  if(image) {
    image(cormat)
    image(as.matrix(ash_cor_only))
  }

  if(all.equal(target=ash_cor_only, current=ash_cor_PD2, tolerance=tol)==TRUE){
    cat("ash cor only and ash cor PD matrices are same")
  }else{
    cat("ash cor only and ash cor PD matrices are different")
  }

  ll <- list("ash_cor_only" = ash_cor_only, "ash_cor_PD" = ash_cor_PD2)
  return(ll)

}

autoselect.mixsd = function(betahat,sebetahat,mult){
  sebetahat=sebetahat[sebetahat!=0] #To avoid exact measure causing (usually by mistake)
  sigmaamin = min(sebetahat)/5 #so that the minimum is small compared with measurement precision
  if(all(betahat^2<=sebetahat^2)){
    sigmaamax = 5*sigmaamin #to deal with the occassional odd case where this could happen; 8 is arbitrary
  }else{
    sigmaamax = 2*sqrt(max(betahat^2-sebetahat^2)) #this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2
  }
  if(mult==0){
    return(c(0,sigmaamax/2))
  }else{
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint):0) * sigmaamax)
  }
}

setprior=function(prior,k,nullweight,null.comp){
  if(!is.numeric(prior)){
    if(prior=="nullbiased"){ # set up prior to favour "null"
      prior = rep(1,k)
      prior[null.comp] = nullweight #prior 10-1 in favour of null by default
    }else if(prior=="uniform"){
      prior = rep(1,k)
    } else if(prior=="unit"){
      prior = rep(1/k,k)
    }
  }
  if(length(prior)!=k | !is.numeric(prior)){
    stop("invalid prior specification")
  }
  return(prior)
}

initpi = function(k,n,null.comp,randomstart){
  if(randomstart){
    pi = stats::rgamma(k,1,1)
  } else {
    if(k<n){
      pi=rep(1,k)/n #default initialization strongly favours null; puts weight 1/n on everything except null
      pi[null.comp] = (n-k+1)/n #the motivation is data can quickly drive away from null, but tend to drive only slowly toward null.
    } else {
      pi=rep(1,k)/k
    }
  }
  pi=normalize(pi)
  return(pi)
}

normalize <- function(x) { return(x/sum(x))}

diriKL = function(p,q){
  p.sum = sum(p)
  q.sum = sum(q)
  k = length(q)
  KL = lgamma(q.sum)-lgamma(p.sum)+sum((q-p)*(digamma(q)-digamma(rep(q.sum,k))))+sum(lgamma(p)-lgamma(q))
  return(KL)
}

estimate_mixprop = function(betahat,sebetahat,g,prior,optmethod=c("mixEM","mixVBEM","cxxMixSquarem","mixIP"),null.comp=1,df=NULL,control=list()){
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  optmethod=match.arg(optmethod)
  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)

  pi_init = g$pi
  if(optmethod=="mixVBEM"){pi_init=NULL}  #for some reason pi_init doesn't work with mixVBEM

  k=ncomp(g)
  n = length(betahat)
  controlinput$tol = min(0.1/n,1.e-7) # set convergence criteria to be more stringent for larger samples

  if(controlinput$trace==TRUE){tic()}

  matrix_llik = t(log_compdens_conv(g,betahat,sebetahat,df)) #an n by k matrix
  matrix_llik = matrix_llik - apply(matrix_llik,1, max) #avoid numerical issues by subtracting max of each row
  matrix_lik = exp(matrix_llik)

  # the last of these conditions checks whether the gradient at the null is negative wrt pi0
  # to avoid running the optimization when the global null (pi0=1) is the optimal.
  if(optmethod=="mixVBEM" || max(prior[-1])>1 || min(gradient(matrix_lik)+prior[1]-1,na.rm=TRUE)<0){
    fit=do.call(optmethod,args = list(matrix_lik= matrix_lik, prior=prior, pi_init=pi_init, control=controlinput))
  } else {
    fit = list(converged=TRUE,pihat=c(1,rep(0,k-1)))
  }

  ## check if IP method returns negative mixing proportions. If so, run EM.
  if (optmethod == "mixIP" & (min(fit$pihat) < -10 ^ -12)) {
    message("Interior point method returned negative mixing proportions.\n Switching to EM optimization.")
    optmethod <- "mixEM"
    fit = do.call(optmethod, args = list(matrix_lik = matrix_lik,
                                         prior = prior, pi_init = pi_init,
                                         control = controlinput))
  }

  if(!fit$converged){
    warning("Optimization failed to converge. Results may be unreliable. Try increasing maxiter and rerunning.")
  }

  pi = fit$pihat
  converged = fit$converged
  niter = fit$niter

  loglik.final =  penloglik(pi,matrix_lik,1) #compute penloglik without penalty
  null.loglik = sum(log(matrix_lik[,null.comp]))
  g$pi=pi
  if(controlinput$trace==TRUE){toc()}

  return(list(loglik=loglik.final,null.loglik=null.loglik,
              matrix_lik=matrix_lik,converged=converged,g=g,niter=niter))
}

mixVBEM = function(matrix_lik, prior, pi_init = NULL,control=list()){
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)

  k=ncol(matrix_lik)
  if(is.null(pi_init)){  pi_init = rep(1,k)  }# Use as starting point for pi
  res = SQUAREM::squarem(par=pi_init,fixptfn=VBfixpoint, objfn=VBnegpenloglik,matrix_lik=matrix_lik, prior=prior, control=controlinput)

  return(list(pihat = res$par/sum(res$par), B=res$value.objfn, niter = res$iter, converged=res$convergence,post=res$par))
}


VBfixpoint = function(pipost, matrix_lik, prior){
  n=nrow(matrix_lik)
  k=ncol(matrix_lik)
  avgpipost = matrix(exp(rep(digamma(pipost),n)-rep(digamma(sum(pipost)),k*n)),ncol=k,byrow=TRUE)
  classprob = avgpipost*matrix_lik
  classprob = classprob/rowSums(classprob) # n by k matrix
  pipostnew = colSums(classprob) + prior
  return(pipostnew)
}

VBnegpenloglik=function(pipost,matrix_lik,prior){
  return(-VBpenloglik(pipost,matrix_lik,prior))
}

VBpenloglik = function(pipost, matrix_lik, prior){
  n=nrow(matrix_lik)
  k=ncol(matrix_lik)
  avgpipost = matrix(exp(rep(digamma(pipost),n)-rep(digamma(sum(pipost)),k*n)),ncol=k,byrow=TRUE)
  classprob = avgpipost*matrix_lik
  classprob = classprob/rowSums(classprob) # n by k matrix

  B= sum(classprob*log(avgpipost*matrix_lik),na.rm=TRUE) - diriKL(prior,pipost) - sum(classprob*log(classprob))
  return(B)
}




compdens = function(m,y,log=FALSE){
  UseMethod("compdens")
}

compdens.default = function(m,y,log=FALSE){
  stop(paste("Invalid class", class(m), "for first argument in",  match.call()))
}

comp_sd = function(m){
  UseMethod("comp_sd")
}
comp_sd.default = function(m){
  stop("method comp_sd not written for this class")
}

comp_mean2 = function(m){
  UseMethod("comp_mean2")
}
comp_mean2.default = function(m){
  comp_sd(m)^2 + comp_mean(m)^2
}


calc_mixmean = function(m){
  UseMethod("calc_mixmean")
}
calc_mixmean.default = function(m){
  sum(m$pi * comp_mean(m))
}


mixmean2 = function(m){
  UseMethod("mixmean2")
}
mixmean2.default = function(m){
  sum(m$pi * comp_mean2(m))
}

calc_mixsd = function(m){
  UseMethod("calc_mixsd")
}
calc_mixsd.default = function(m){
  sqrt(mixmean2(m)-calc_mixmean(m)^2)
}

comp_mean = function(m){
  UseMethod("comp_mean")
}
comp_mean.default = function(m){
  stop("method comp_mean not written for this class")
}


ncomp = function(m){
  UseMethod("ncomp")
}

ncomp.default = function(m){
  return(length(m$pi))
}


mixprop = function(m){
  UseMethod("mixprop")
}
mixprop.default = function(m){
  m$pi
}


mixcdf = function(x,y,lower.tail=TRUE){
  UseMethod("mixcdf")
}

mixcdf.default = function(x,y,lower.tail=TRUE){
  x$pi %*% comp_cdf(x,y,lower.tail)
}

comp_cdf = function(m,y,lower.tail=TRUE){
  UseMethod("comp_cdf")
}

comp_cdf.default = function(m,y,lower.tail=TRUE){
  stop("comp_cdf not implemented for this class")
}


dens = function(x,y){
  UseMethod("dens")
}

dens.default = function(x,y){
  return (x$pi %*% compdens(x, y))
}

#find log likelihood of data in x (a vector) for mixture in m
loglik = function(m,x){
  UseMethod("loglik")
}
loglik.default = function(m,x){
  sum(log(dens(m,x)))
}

loglik_conv = function(m,betahat,betahatsd,v,FUN="+"){
  UseMethod("loglik_conv")
}

loglik_conv.default = function(m,betahat,betahatsd,v,FUN="+"){
  sum(log(dens_conv(m,betahat,betahatsd,v,FUN)))
}


compdens_conv = function(m,x,s,v,FUN="+"){
  UseMethod("compdens_conv")
}
compdens_conv.default = function(m,x,s,v,FUN="+"){
  stop(paste("Invalid class", class(m), "for first argument in",  match.call()))
}

log_compdens_conv = function(m,x,s,v,FUN="+"){
  UseMethod("log_compdens_conv")
}
log_compdens_conv.default = function(m,x,s,v,FUN="+"){
  log(compdens_conv(m,x,s,v,FUN))
}



dens_conv = function(m,x,s,v,FUN="+"){
  UseMethod("dens_conv")
}
dens_conv.default = function(m,x,s,v,FUN="+"){
  colSums(m$pi * compdens_conv(m,x,s,v,FUN))
}


comppostprob=function(m,x,s,v){
  UseMethod("comppostprob")
}


old.comppostprob.default = function(m,x,s,v){
  tmp= (t(m$pi * compdens_conv(m,x,s,v))/dens_conv(m,x,s,v))
  ismissing = (is.na(x) | is.na(s))
  tmp[ismissing,]=m$pi
  t(tmp)
}

comppostprob.default = function(m,x,s,v){
  lpost = log_compdens_conv(m,x,s,v) + log(m$pi) # lpost is k by n of log posterior prob (unnormalized)
  lpmax = apply(lpost,2,max) #dmax is of length n
  tmp = exp(t(lpost)-lpmax) #subtracting the max of the logs is just done for numerical stability
  tmp = tmp/rowSums(tmp)
  ismissing = (is.na(x) | is.na(s))
  tmp[ismissing,]=m$pi
  t(tmp)
}


compcdf_post=function(m,c,betahat,sebetahat,v){
  UseMethod("compcdf_post")
}

compcdf_post.default=function(m,c,betahat,sebetahat,v){
  stop("method compcdf_post not written for this class")
}

cdf_post = function(m,c,betahat,sebetahat,v){
  UseMethod("cdf_post")
}

cdf_post.default=function(m,c,betahat,sebetahat,v){
  colSums(comppostprob(m,betahat,sebetahat,v)*compcdf_post(m,c,betahat,sebetahat,v))
}

vcdf_post = function(m,c,betahat,sebetahat,v){
  mapply(cdf_post,c,MoreArgs=list(m=m,betahat=betahat,sebetahat=sebetahat,v=v))
}

pcdf_post = function(m,c,betahat,sebetahat,v){
  mapply(cdf_post,c,betahat,sebetahat,MoreArgs=list(m=m,v=v))
}

postmean = function(m, betahat,sebetahat,v){
  UseMethod("postmean")
}


postmean.default = function(m,betahat,sebetahat,v){
  colSums(comppostprob(m,betahat,sebetahat,v) * comp_postmean(m,betahat,sebetahat,v))
}


postmean2 = function(m, betahat,sebetahat,v){
  UseMethod("postmean2")
}

postmean2.default = function(m,betahat,sebetahat,v){
  colSums(comppostprob(m,betahat,sebetahat,v) * comp_postmean2(m,betahat,sebetahat,v))
}

postsd = function(m,betahat,sebetahat,v){
  UseMethod("postsd")
}


postsd.default = function(m,betahat,sebetahat,v){
  sqrt(postmean2(m,betahat,sebetahat,v)-postmean(m,betahat,sebetahat,v)^2)
}

comp_postmean2 = function(m,betahat,sebetahat,v){
  UseMethod("comp_postmean2")
}
comp_postmean2.default = function(m,betahat,sebetahat,v){
  stop("method comp_postmean2 not written for this class")
  #comp_postsd(m,betahat,sebetahat,v)^2 + comp_postmean(m,betahat,sebetahat,v)^2
}

comp_postmean = function(m, betahat,sebetahat,v){
  UseMethod("comp_postmean")
}
comp_postmean.default = function(m,betahat,sebetahat,v){
  stop("method comp_postmean not written for this class")
}


comp_postsd = function(m, betahat,sebetahat,v){
  UseMethod("comp_postsd")
}
comp_postsd.default = function(m,betahat,sebetahat,v){
  stop("method comp_postsd not written for this class")
}

min_lim = function(m){
  UseMethod("min_lim")
}
min_lim.default=function(m){
  -5
}

max_lim = function(m){
  UseMethod("max_lim")
}
max_lim.default=function(m){
  5
}


#plot density of mixture
plot_dens = function(m,npoints=100,...){
  UseMethod("plot_dens")
}
plot_dens.default = function(m,npoints=100,...){
  x = seq(min_lim(m),max_lim(m),length=npoints)
  graphics::plot(x,dens(m,x),type="l",xlab="density",ylab="x",...)
}

plot_post_cdf = function(m,betahat,sebetahat,v,npoints=100,...){
  UseMethod("plot_post_cdf")
}
plot_post_cdf.default = function(m,betahat,sebetahat,v,npoints=100,...){
  x = seq(min_lim(m),max_lim(m),length=npoints)
  x_cdf = vapply(x,cdf_post,FUN.VALUE=betahat,m=m,betahat=betahat,sebetahat=sebetahat,v=v)
  graphics::plot(x,x_cdf,type="l",xlab="x",ylab="cdf",...)
  # for(i in 2:nrow(x_cdf)){
  #   lines(x,x_cdf[i,],col=i)
  # }
}


normalmix = function(pi,mean,sd){
  structure(data.frame(pi,mean,sd),class="normalmix")
}


comp_sd.normalmix = function(m){
  m$sd
}

comp_mean.normalmix = function(m){
  m$mean
}

compdens.normalmix = function(m,y,log=FALSE){
  k=ncomp(m)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(stats::dnorm(d, m$mean, m$sd, log),nrow=k))
}

compdens_conv.normalmix = function(m,x,s,v,FUN="+"){
  if(!is.null(v)){
    stop("method compdens_conv of normal mixture not written for df!=NULL")
  }
  if(length(s)==1){s=rep(s,length(x))}
  sdmat = sqrt(outer(s^2,m$sd^2,FUN)) #n by k matrix of standard deviations of convolutions
  return(t(stats::dnorm(outer(x,m$mean,FUN="-")/sdmat)/sdmat))
}


log_compdens_conv.normalmix = function(m,x,s,v,FUN="+"){
  if(!is.null(v)){
    stop("method compdens_conv of normal mixture not written for df!=NULL")
  }
  if(length(s)==1){s=rep(s,length(x))}
  sdmat = sqrt(outer(s^2,m$sd^2,FUN)) #n by k matrix of standard deviations of convolutions
  return(t(stats::dnorm(outer(x,m$mean,FUN="-")/sdmat,log=TRUE) - log(sdmat)))
}


comp_cdf.normalmix = function(m,y,lower.tail=TRUE){
  vapply(y,stats::pnorm,m$mean,m$mean,m$sd,lower.tail)
}

compcdf_post.normalmix=function(m,c,betahat,sebetahat,v){
  if(!is.null(v)){
    stop("Error: normal mixture for student-t likelihood is not yet implemented")
  }
  k = length(m$pi)
  n=length(betahat)
  #compute posterior standard deviation (s1) and posterior mean (m1)
  s1 = sqrt(outer(sebetahat^2,m$sd^2,FUN="*")/outer(sebetahat^2,m$sd^2,FUN="+"))
  ismissing = (is.na(betahat) | is.na(sebetahat))
  s1[ismissing,]=m$sd

  m1 = t(comp_postmean(m,betahat,sebetahat,v))
  t(stats::pnorm(c,mean=m1,sd=s1))
}


comp_postmean.normalmix = function(m,betahat,sebetahat,v){
  if(!is.null(v)){
    stop("method comp_postmean of normal mixture not written for df!=NULL")
  }
  tmp=(outer(sebetahat^2,m$mean, FUN="*") + outer(betahat,m$sd^2, FUN="*"))/outer(sebetahat^2,m$sd^2,FUN="+")
  ismissing = (is.na(betahat) | is.na(sebetahat))
  tmp[ismissing,]=m$mean #return prior mean when missing data
  t(tmp)
}

comp_postsd.normalmix = function(m,betahat,sebetahat,v){
  if(!is.null(v)){
    stop("method comp_postsd of normal mixture not written for df!=NULL")
  }
  t(sqrt(outer(sebetahat^2,m$sd^2,FUN="*")/outer(sebetahat^2,m$sd^2,FUN="+")))
}

comp_postmean2.normalmix = function(m,betahat,sebetahat,v){
  comp_postsd(m,betahat,sebetahat,v)^2 + comp_postmean(m,betahat,sebetahat,v)^2
}

loglik_conv_mixlik = function(m,betahat,betahatsd,v,pilik,FUN="+"){
  UseMethod("loglik_conv_mixlik")
}

loglik_conv_mixlik.default = function(m,betahat,betahatsd,v,pilik,FUN="+"){
  sum(log(dens_conv_mixlik(m,betahat,betahatsd,v,pilik,FUN)))
}


compdens_conv_mixlik = function(m,x,s,v,pilik,FUN="+"){
  UseMethod("compdens_conv_mixlik")
}
compdens_conv_mixlik.default = function(m,x,s,v,pilik,FUN="+"){
  dens=NULL
  for (i in 1:dim(pilik)[2]){
    dens=rbind(dens,pilik[,i]*compdens_conv(m,x,s[,i],v[i],FUN))
  }
  return(dens)
}


dens_conv_mixlik = function(m,x,s,v,pilik,FUN="+"){
  UseMethod("dens_conv_mixlik")
}
dens_conv_mixlik.default = function(m,x,s,v,pilik,FUN="+"){
  l=dim(pilik)[2]
  colSums(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik,FUN))
}


comppostprob_mixlik=function(m,x,s,v,pilik){
  UseMethod("comppostprob_mixlik")
}
comppostprob_mixlik.default = function(m,x,s,v,pilik){
  l=dim(pilik)[2]
  k=length(m$pi)
  tmp= (t(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik))/dens_conv_mixlik(m,x,s,v,pilik))
  ismissing = (is.na(x) | apply(is.na(s),1,sum))
  tmp[ismissing,]=rep(m$pi,l)/l
  group=rep(1:k,l)
  return(rowsum(t(tmp),group))
}



comppostprob_mixlik2=function(m,x,s,v,pilik){
  UseMethod("comppostprob_mixlik2")
}

comppostprob_mixlik2.default = function(m,x,s,v,pilik){
  l=dim(pilik)[2]
  k=length(m$pi)
  tmp= (t(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik))/dens_conv_mixlik(m,x,s,v,pilik))
  ismissing = (is.na(x) | apply(is.na(s),1,sum))
  tmp[ismissing,]=rep(m$pi,l)/l
  return(t(tmp))
}



compcdf_post_mixlik=function(m,c,betahat,sebetahat,v,pilik){
  UseMethod("compcdf_post_mixlik")
}

compcdf_post_mixlik.default=function(m,c,betahat,sebetahat,v,pilik){
  cdf=NULL
  for (i in 1:dim(pilik)[2]){
    cdf=rbind(cdf,pilik[,i]*compcdf_post(m,c,betahat,sebetahat[,i],v[i]))
  }
  cdf
}

cdf_post_mixlik = function(m,c,betahat,sebetahat,v,pilik){
  UseMethod("cdf_post_mixlik")
}
cdf_post_mixlik.default=function(m,c,betahat,sebetahat,v,pilik){
  colSums(comppostprob_mixlik2(m,betahat,sebetahat,v,pilik)*
            compcdf_post_mixlik(m,c,betahat,sebetahat,v,pilik))
}

postmean_mixlik = function(m, betahat,sebetahat,v,pilik){
  UseMethod("postmean_mixlik")
}
postmean_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  colSums(comppostprob_mixlik2(m,betahat,sebetahat,v,pilik) * comp_postmean_mixlik(m,betahat,sebetahat,v,pilik))
}


postmean2_mixlik = function(m, betahat,sebetahat,v,pilik){
  UseMethod("postmean2_mixlik")
}

postmean2_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  colSums(comppostprob_mixlik2(m,betahat,sebetahat,v,pilik) * comp_postmean2_mixlik(m,betahat,sebetahat,v,pilik))
}

postsd_mixlik = function(m,betahat,sebetahat,v,pilik){
  UseMethod("postsd_mixlik")
}
postsd_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  sqrt(postmean2_mixlik(m,betahat,sebetahat,v,pilik)-postmean_mixlik(m,betahat,sebetahat,v,pilik)^2)
}


comp_postmean2_mixlik = function(m,betahat,sebetahat,v,pilik){
  UseMethod("comp_postmean2_mixlik")
}

comp_postmean2_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  comp_postsd_mixlik(m,betahat,sebetahat,v,pilik)^2 +
    comp_postmean_mixlik(m,betahat,sebetahat,v,pilik)^2
}


comp_postmean_mixlik=function(m,betahat,sebetahat,v,pilik){
  UseMethod("comp_postmean_mixlik")
}
comp_postmean_mixlik.default=function(m,betahat,sebetahat,v,pilik){
  mean=NULL
  for (i in 1:dim(pilik)[2]){
    mean=rbind(mean,comp_postmean(m,betahat,sebetahat[,i],v[i]))
  }
  return(mean)
}



comp_postsd_mixlik=function(m,betahat,sebetahat,v,pilik){
  UseMethod("comp_postsd_mixlik")
}
comp_postsd_mixlik.default=function(m,betahat,sebetahat,v,pilik){
  sd=NULL
  for (i in 1:dim(pilik)[2]){
    sd=rbind(sd,comp_postsd(m,betahat,sebetahat[,i],v[i]))
  }
  return(sd)
}

negpenloglik = function(pi,matrix_lik,prior){return(-penloglik(pi,matrix_lik,prior))}

penloglik = function(pi, matrix_lik, prior){
  pi = normalize(pmax(0,pi))
  m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  loglik = sum(log(m.rowsum))
  subset = (prior != 1.0)
  priordens = sum((prior-1)[subset]*log(pi[subset]))
  return(loglik+priordens)
}

