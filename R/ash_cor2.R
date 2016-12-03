#' @title Adaptive shrinkage of a correlation matrix with Variational EM
#'
#' @description This function performs adaptive shrinkage of a correlation matrix
#'
#' @param cormat The estimated sample correlation matrix
#' @param nsamples The number of samples over which the correlation matrix is estimated.
#' @param image if TRUE, plots an image of the shrunk and non-shrunk correlation matrices
#' @param tol The tolerance to check the difference between ash-cor only and ash-cor PD matrices.
#' @param nullweight The weight of the null component of the ash prior.
#' @param null.comp The number of null components. Default is 1.
#'
#' @return Returns a shrunk version of the sample correlation matrix (the output is also a correlation matrix)
#'
#' @keywords adaptive shrinkage, correlation
#' @export

ash_cor2 <- function(cormat, nsamples, image=FALSE, tol=1e-06,
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



################################## GENERIC FUNCTIONS ############################

#' Generic function of calculating the component densities of the
#'     mixture
#' @param m mixture of k components generated by normalmix() or
#'     unimix() or igmix()
#' @param y is an n-vector of location
#' @param log whether to use log-scale on densities
#' @return A k by n matrix of densities
#' @export
compdens = function(m,y,log=FALSE){
  UseMethod("compdens")
}
#' @export
compdens.default = function(m,y,log=FALSE){
  stop(paste("Invalid class", class(m), "for first argument in",  match.call()))
}

#' Generic function of extracting the standard deviations of
#'     components of the mixture
#' @param m a mixture of k components generated by normalmix() or
#'     unimix() or igmix()
#' @return it returns a vector of standard deviations
#' @export
comp_sd = function(m){
  UseMethod("comp_sd")
}
#' @export
comp_sd.default = function(m){
  stop("method comp_sd not written for this class")
}

#' Generic function of calculating the second moment of components of
#'     the mixture
#' @param m a mixture of k components generated by normalmix() or
#'     unimix() or igmix()
#' @return it returns a vector of second moments.
comp_mean2 = function(m){
  UseMethod("comp_mean2")
}
comp_mean2.default = function(m){
  comp_sd(m)^2 + comp_mean(m)^2
}


#' Generic function of calculating the overall mean of the mixture
#'
#'
#' @param m a mixture of k components generated by normalmix() or
#'     unimix() or igmix()
#' @return it returns scalar, the mean of the mixture distribution.
calc_mixmean = function(m){
  UseMethod("calc_mixmean")
}
calc_mixmean.default = function(m){
  sum(m$pi * comp_mean(m))
}


#' Generic function of calculating the overall second moment of the
#'     mixture
#' @param m a mixture of k components generated by normalmix() or
#'     unimix() or igmix()
#' @return it returns scalar
mixmean2 = function(m){
  UseMethod("mixmean2")
}
mixmean2.default = function(m){
  sum(m$pi * comp_mean2(m))
}

#' Generic function of calculating the overall standard deviation of
#'     the mixture
#' @param m a mixture of k components generated by normalmix() or
#'     unimix() or igmix()
#' @return it returns scalar
#' @export
calc_mixsd = function(m){
  UseMethod("calc_mixsd")
}
calc_mixsd.default = function(m){
  sqrt(mixmean2(m)-calc_mixmean(m)^2)
}

#' Generic function of calculating the first moment of components of
#'     the mixture
#' @param m a mixture of k components generated by normalmix() or
#'     unimix() or igmix()
#' @return it returns a vector of means.
#' @export
comp_mean = function(m){
  UseMethod("comp_mean")
}
#' @export
comp_mean.default = function(m){
  stop("method comp_mean not written for this class")
}

#number of components
#' ncomp
#' @param m a mixture of k components generated by normalmix() or
#'     unimix() or igmix()
#' @export
#'
ncomp = function(m){
  UseMethod("ncomp")
}

#' @title ncomp.default
#' @description The default version of \code{\link{ncomp}}.
#' @param m a mixture of k components generated by normalmix() or
#'     unimix() or igmix()
#' @export
#'
ncomp.default = function(m){
  return(length(m$pi))
}


#' Generic function of extracting the mixture proportions
#' @param m a mixture of k components generated by normalmix() or
#'     unimix() or igmix()
#' @return it returns a vector of component probabilities, summing up
#'     to 1.
#' @export
mixprop = function(m){
  UseMethod("mixprop")
}
#' @export
mixprop.default = function(m){
  m$pi
}

#' @title mixcdf
#'
#' @description Returns cdf for a mixture (generic function)
#'
#' @details None
#'
#' @param x a mixture (eg of type normalmix or unimix)
#' @param y locations at which cdf to be computed
#' @param lower.tail boolean indicating whether to report lower tail
#'
#' @return an object of class normalmix
#'
#' @export
#'
#' @examples mixcdf(normalmix(c(0.5,0.5),c(0,0),c(1,2)),seq(-4,4,length=100))
#'
mixcdf = function(x,y,lower.tail=TRUE){
  UseMethod("mixcdf")
}
#' @title mixcdf.default
#' @description The default version of \code{\link{mixcdf}}.
#' @param x a mixture (eg of type normalmix or unimix)
#' @param y locations at which cdf to be computed
#' @param lower.tail boolean indicating whether to report lower tail
#' @export
#'
mixcdf.default = function(x,y,lower.tail=TRUE){
  x$pi %*% comp_cdf(x,y,lower.tail)
}

#find cdf for each component, a generic function
#' Generic function of computing the cdf for each component
#' @param m a mixture (eg of type normalmix or unimix)
#' @param y locations at which cdf to be computed
#' @param lower.tail boolean indicating whether to report lower tail
#' @return it returns a vector of probabilities, with length equals to
#'     number of components in m
#' @export
comp_cdf = function(m,y,lower.tail=TRUE){
  UseMethod("comp_cdf")
}
#' @export
comp_cdf.default = function(m,y,lower.tail=TRUE){
  stop("comp_cdf not implemented for this class")
}


#' Find density at y, a generic function
#' @param x A mixture of k components generated by
#'     \code{\link{normalmix}} or \code{\link{unimix}}.
#' @param y An n-vector of the location.
#' @export
dens = function(x,y){
  UseMethod("dens")
}
#' @export
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


#' @title loglik_conv
#' @description find log likelihood of data in betahat, when the
#'     mixture m is convolved with a normal with sd betahatsd
#'
#' @param m mixture distribution with k components
#' @param betahat an n vector of observations
#' @param betahatsd an n vector of standard errors
#' @param v degree of freedom of error distribution
#' @param FUN default is "+"
#' @export
#'
loglik_conv = function(m,betahat,betahatsd,v,FUN="+"){
  UseMethod("loglik_conv")
}
#' @title loglik_conv.default
#'
#' @description The default version of \code{\link{loglik_conv}}.
#'
#' @param m mixture distribution with k components
#' @param betahat an n vector of observations
#' @param betahatsd an n vector of standard errors
#' @param v degree of freedom of error distribution
#' @param FUN default is "+"
#' @export
#'
loglik_conv.default = function(m,betahat,betahatsd,v,FUN="+"){
  sum(log(dens_conv(m,betahat,betahatsd,v,FUN)))
}



#' @title compdens_conv
#' @description compute the density of the components of the mixture m
#'     when convoluted with a normal with standard deviation s or a
#'     scaled (se) student.t with df v, the density is evaluated at x
#' @param m mixture distribution
#' @param x an n vector
#' @param s an n vector or integer
#' @param v degree of freedom of error distribution
#' @param FUN default is "+"
#' @return a k by n matrix of densities
compdens_conv = function(m,x,s,v,FUN="+"){
  UseMethod("compdens_conv")
}
compdens_conv.default = function(m,x,s,v,FUN="+"){
  stop(paste("Invalid class", class(m), "for first argument in",  match.call()))
}

#' @title log_compdens_conv
#' @description compute the log density of the components of the
#'     mixture m when convoluted with a normal with standard deviation
#'     s or a scaled (se) student.t with df v, the density is
#'     evaluated at x
#' @param m mixture distribution
#' @param x an n vector
#' @param s an n vector or integer
#' @param v degree of freedom of error distribution
#' @param FUN default is "+"
#' @return a k by n matrix of log densities
log_compdens_conv = function(m,x,s,v,FUN="+"){
  UseMethod("log_compdens_conv")
}
log_compdens_conv.default = function(m,x,s,v,FUN="+"){
  log(compdens_conv(m,x,s,v,FUN))
}



#' @title dens_conv
#' @description compute density of mixture m convoluted with normal of
#'     sd (s) or student t with df v at locations x
#' @param m mixture distribution
#' @param x an n vector
#' @param s an n vector or integer
#' @param v degree of freedom of error distribution
#' @param FUN default is "+"
dens_conv = function(m,x,s,v,FUN="+"){
  UseMethod("dens_conv")
}
dens_conv.default = function(m,x,s,v,FUN="+"){
  colSums(m$pi * compdens_conv(m,x,s,v,FUN))
}


#' @title comppostprob
#' @description compute the posterior prob that each observation came
#'     from each component of the mixture m,output a k by n vector of
#'     probabilities computed by weighting the component densities by
#'     pi and then normalizing
#'
#' @param m mixture distribution
#' @param x an n vector
#' @param s an n vector or integer
#' @param v degree of freedom of error distribution
#' @export
comppostprob=function(m,x,s,v){
  UseMethod("comppostprob")
}

#' The old default version of \code{\link{comppostprob}}.
#'
#' @inheritParams comppostprob
#'
#' @export
#'
old.comppostprob.default = function(m,x,s,v){
  tmp= (t(m$pi * compdens_conv(m,x,s,v))/dens_conv(m,x,s,v))
  ismissing = (is.na(x) | is.na(s))
  tmp[ismissing,]=m$pi
  t(tmp)
}

#' @export
comppostprob.default = function(m,x,s,v){
  lpost = log_compdens_conv(m,x,s,v) + log(m$pi) # lpost is k by n of log posterior prob (unnormalized)
  lpmax = apply(lpost,2,max) #dmax is of length n
  tmp = exp(t(lpost)-lpmax) #subtracting the max of the logs is just done for numerical stability
  tmp = tmp/rowSums(tmp)
  ismissing = (is.na(x) | is.na(s))
  tmp[ismissing,]=m$pi
  t(tmp)
}

#' @title compcdf_post
#' @description evaluate cdf of posterior distribution of beta at c. m
#'     is the prior on beta, a mixture; c is location of evaluation
#'     assumption is betahat | beta ~ t_v(beta,sebetahat)
#' @param m mixture distribution with k components
#' @param c a scalar
#' @param betahat an n vector of observations
#' @param sebetahat an n vector of standard errors
#' @param v degree of freedom of error distribution
#' @return a k by n matrix
#' @examples
#' beta = rnorm(100,0,1)
#' betahat= beta+rnorm(100,0,1)
#' sebetahat=rep(1,100)
#' ash.beta = ash(betahat,1,mixcompdist="normal")
#' compcdf_post(ash.beta$fitted.g,0,betahat,sebetahat,NULL)
#' @export
compcdf_post=function(m,c,betahat,sebetahat,v){
  UseMethod("compcdf_post")
}
#' @export
compcdf_post.default=function(m,c,betahat,sebetahat,v){
  stop("method compcdf_post not written for this class")
}

#' @title cdf_post
#' @description evaluate cdf of posterior distribution of beta at c. m
#'     is the prior on beta, a mixture; c is location of evaluation
#'     assumption is betahat | beta ~ t_v(beta,sebetahat)
#' @param m mixture distribution with k components
#' @param c a scalar
#' @param betahat an n vector of observations
#' @param sebetahat an n vector of standard errors
#' @param v degree of freedom of error distribution
#' @return an n vector containing the cdf for beta_i at c
#' @examples
#' beta = rnorm(100,0,1)
#' betahat= beta+rnorm(100,0,1)
#' sebetahat=rep(1,100)
#' ash.beta = ash(betahat,1,mixcompdist="normal")
#' cdf0 = cdf_post(ash.beta$fitted.g,0,betahat,sebetahat,NULL)
#' graphics::plot(cdf0,1-ash.beta$PositiveProb)
#' @export
cdf_post = function(m,c,betahat,sebetahat,v){
  UseMethod("cdf_post")
}
#' @export
cdf_post.default=function(m,c,betahat,sebetahat,v){
  colSums(comppostprob(m,betahat,sebetahat,v)*compcdf_post(m,c,betahat,sebetahat,v))
}
#' @title vcdf_post
#' @description vectorized version of \code{\link{cdf_post}}
#' @param m mixture distribution with k components
#' @param c a numeric vector
#' @param betahat an n vector of observations
#' @param sebetahat an n vector of standard errors
#' @param v degree of freedom of error distribution
#' @return an n vector containing the cdf for beta_i at c
#' @examples
#' beta = rnorm(100,0,1)
#' betahat= beta+rnorm(100,0,1)
#' sebetahat=rep(1,100)
#' ash.beta = ash(betahat,1,mixcompdist="normal")
#' c = vcdf_post(ash.beta$fitted.g,seq(-5,5,length=1000),betahat,sebetahat,NULL)
#' @export
vcdf_post = function(m,c,betahat,sebetahat,v){
  mapply(cdf_post,c,MoreArgs=list(m=m,betahat=betahat,sebetahat=sebetahat,v=v))
}
#' @title pcdf_post
#' @description ``parallel" vector version of \code{\link{cdf_post}} where c is a vector, of same length as betahat and sebetahat
#' @param m mixture distribution with k components
#' @param c a numeric vector with n elements
#' @param betahat an n vector of observations
#' @param sebetahat an n vector of standard errors
#' @param v degree of freedom of error distribution (scalar)
#' @return an n vector, whose ith element is the cdf for beta_i at c_i
#' @examples
#' beta = rnorm(100,0,1)
#' betahat= beta+rnorm(100,0,1)
#' sebetahat=rep(1,100)
#' ash.beta = ash(betahat,1,mixcompdist="normal")
#' c = pcdf_post(ash.beta$fitted.g,beta,betahat,sebetahat,NULL)
#' @export
pcdf_post = function(m,c,betahat,sebetahat,v){
  mapply(cdf_post,c,betahat,sebetahat,MoreArgs=list(m=m,v=v))
}

#output posterior mean for beta for prior mixture m,
#given observations betahat, sebetahat, df v
#' postmean
#'
#' @param m mixture distribution with k components
#' @param betahat an n vector of observations
#' @param sebetahat an n vector of standard errors
#' @param v degree of freedom of error distribution
#' @export
postmean = function(m, betahat,sebetahat,v){
  UseMethod("postmean")
}

#' @export
postmean.default = function(m,betahat,sebetahat,v){
  colSums(comppostprob(m,betahat,sebetahat,v) * comp_postmean(m,betahat,sebetahat,v))
}



#' @title postmean2
#' @description output posterior mean-squared value for beta for prior
#'     mixture m,given observations betahat, sebetahat, df v
#' @param m mixture distribution with k components
#' @param betahat an n vector of observations
#' @param sebetahat an n vector of standard errors
#' @param v degree of freedom of error distribution
#' @export
postmean2 = function(m, betahat,sebetahat,v){
  UseMethod("postmean2")
}
#' @export
postmean2.default = function(m,betahat,sebetahat,v){
  colSums(comppostprob(m,betahat,sebetahat,v) * comp_postmean2(m,betahat,sebetahat,v))
}


#' @title postsd
#' @description output posterior sd for beta for prior mixture m,given
#'     observations betahat, sebetahat, df v
#' @param m mixture distribution with k components
#' @param betahat an n vector of observations
#' @param sebetahat an n vector of standard errors
#' @param v degree of freedom of error distribution
#' @export
postsd = function(m,betahat,sebetahat,v){
  UseMethod("postsd")
}

#' @export
postsd.default = function(m,betahat,sebetahat,v){
  sqrt(postmean2(m,betahat,sebetahat,v)-postmean(m,betahat,sebetahat,v)^2)
}



#' @title comp_postmean2
#' @description output posterior mean-squared value for beta for prior
#'     mixture m,given observations betahat, sebetahat, df v
#' @param m mixture distribution with k components
#' @param betahat an n vector of observations
#' @param sebetahat an n vector of standard errors
#' @param v degree of freedom of error distribution
#' @export
comp_postmean2 = function(m,betahat,sebetahat,v){
  UseMethod("comp_postmean2")
}
#' @export
comp_postmean2.default = function(m,betahat,sebetahat,v){
  stop("method comp_postmean2 not written for this class")
  #comp_postsd(m,betahat,sebetahat,v)^2 + comp_postmean(m,betahat,sebetahat,v)^2
}



#' @title comp_postmean
#' @description output posterior mean for beta for each component of
#'     prior mixture m,given observations betahat, sebetahat, df v
#' @param m mixture distribution with k components
#' @param betahat an n vector of observations
#' @param sebetahat an n vector of standard errors
#' @param v degree of freedom of error distribution
#' @export
comp_postmean = function(m, betahat,sebetahat,v){
  UseMethod("comp_postmean")
}
#' @export
comp_postmean.default = function(m,betahat,sebetahat,v){
  stop("method comp_postmean not written for this class")
}



#' @title comp_postsd
#' @description output posterior sd for beta for each component of
#'     prior mixture m,given observations betahat, sebetahat, df v
#' @param m mixture distribution with k components
#' @param betahat an n vector of observations
#' @param sebetahat an n vector of standard errors
#' @param v degree of freedom of error distribution
#' @examples
#' beta = rnorm(100,0,1)
#' betahat= beta+rnorm(100,0,1)
#' ash.beta = ash(betahat,1,mixcompdist="normal")
#' comp_postmean(ash.beta$fitted.g,betahat,rep(1,100),NULL)
#' comp_postsd(ash.beta$fitted.g,betahat,rep(1,100),NULL)
#' comppostprob(ash.beta$fitted.g,betahat,rep(1,100),NULL)
#' @export
comp_postsd = function(m, betahat,sebetahat,v){
  UseMethod("comp_postsd")
}
#' @export
comp_postsd.default = function(m,betahat,sebetahat,v){
  stop("method comp_postsd not written for this class")
}

#find nice limits of mixture m for plotting
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

############################### METHODS FOR normalmix class ###########################

#' @title Constructor for normalmix class
#'
#' @description Creates an object of class normalmix (finite mixture
#'     of univariate normals)
#'
#' @details None
#'
#' @param pi vector of mixture proportions
#' @param mean vector of means
#' @param sd vector of standard deviations
#'
#' @return an object of class normalmix
#'
#' @export
#'
#' @examples normalmix(c(0.5,0.5),c(0,0),c(1,2))
#'
normalmix = function(pi,mean,sd){
  structure(data.frame(pi,mean,sd),class="normalmix")
}


#' @title comp_sd.normalmix
#' @description returns sds of the normal mixture
#' @param m a normal mixture distribution with k components
#' @return a vector of length k
#' @export
comp_sd.normalmix = function(m){
  m$sd
}

#' @title comp_mean.normalmix
#' @description returns mean of the normal mixture
#' @param m a normal mixture distribution with k components
#' @return a vector of length k
#' @export
comp_mean.normalmix = function(m){
  m$mean
}

#' @export
compdens.normalmix = function(m,y,log=FALSE){
  k=ncomp(m)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(stats::dnorm(d, m$mean, m$sd, log),nrow=k))
}


#' @title compdens_conv.normalmix
#' @description returns density of convolution of each component of a
#'     normal mixture with N(0,s^2) or s*t(v) at x. Note that
#'     convolution of two normals is normal, so it works that way
#' @param m mixture distribution with k components
#' @param x an n-vector at which density is to be evaluated
#' @param s an n vector of standard errors
#' @param v degree of freedom of error distribution
#' @param FUN default is "+"
#' @return a k by n matrix
compdens_conv.normalmix = function(m,x,s,v,FUN="+"){
  if(!is.null(v)){
    stop("method compdens_conv of normal mixture not written for df!=NULL")
  }
  if(length(s)==1){s=rep(s,length(x))}
  sdmat = sqrt(outer(s^2,m$sd^2,FUN)) #n by k matrix of standard deviations of convolutions
  return(t(stats::dnorm(outer(x,m$mean,FUN="-")/sdmat)/sdmat))
}

#' @title log_compdens_conv.normalmix
#' @description returns log-density of convolution of each component
#'     of a normal mixture with N(0,s^2) or s*t(v) at x. Note that
#'     convolution of two normals is normal, so it works that way
#' @param m mixture distribution with k components
#' @param x an n-vector at which density is to be evaluated
#' @param s an n vector of standard errors
#' @param v degree of freedom of error distribution
#' @param FUN default is "+"
#' @return a k by n matrix
log_compdens_conv.normalmix = function(m,x,s,v,FUN="+"){
  if(!is.null(v)){
    stop("method compdens_conv of normal mixture not written for df!=NULL")
  }
  if(length(s)==1){s=rep(s,length(x))}
  sdmat = sqrt(outer(s^2,m$sd^2,FUN)) #n by k matrix of standard deviations of convolutions
  return(t(stats::dnorm(outer(x,m$mean,FUN="-")/sdmat,log=TRUE) - log(sdmat)))
}



#' @export
comp_cdf.normalmix = function(m,y,lower.tail=TRUE){
  vapply(y,stats::pnorm,m$mean,m$mean,m$sd,lower.tail)
}




#' @export
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


#' @export
comp_postmean.normalmix = function(m,betahat,sebetahat,v){
  if(!is.null(v)){
    stop("method comp_postmean of normal mixture not written for df!=NULL")
  }
  tmp=(outer(sebetahat^2,m$mean, FUN="*") + outer(betahat,m$sd^2, FUN="*"))/outer(sebetahat^2,m$sd^2,FUN="+")
  ismissing = (is.na(betahat) | is.na(sebetahat))
  tmp[ismissing,]=m$mean #return prior mean when missing data
  t(tmp)
}


#' @export
comp_postsd.normalmix = function(m,betahat,sebetahat,v){
  if(!is.null(v)){
    stop("method comp_postsd of normal mixture not written for df!=NULL")
  }
  t(sqrt(outer(sebetahat^2,m$sd^2,FUN="*")/outer(sebetahat^2,m$sd^2,FUN="+")))
}

#' @export
comp_postmean2.normalmix = function(m,betahat,sebetahat,v){
  comp_postsd(m,betahat,sebetahat,v)^2 + comp_postmean(m,betahat,sebetahat,v)^2
}

#' @title loglik_conv_mixlik
#'
#' @description find log likelihood of data, when the mixture m is convolved with l-comp normal mixture in betahat with mixture sd betahatsd, mixture proportion pilik.
#'
#' @param m mixture distribution
#' @param betahat an n vector
#' @param betahatsd an n-by-l matrix
#' @param v degree of freedom of error distribution
#' @param pilik an l vector
#' @param FUN default is "+"
#' @export
#'
loglik_conv_mixlik = function(m,betahat,betahatsd,v,pilik,FUN="+"){
  UseMethod("loglik_conv_mixlik")
}
#' @title loglik_conv.default
#'
#' @description The default version of \code{\link{loglik_conv}}.
#'
#' @param m mixture distribution
#' @param betahat an n vector
#' @param betahatsd an n-by-l matrix
#' @param v degree of freedom of error distribution
#' @param pilik an l vector
#' @param FUN default is "+"
#' @export
#'
loglik_conv_mixlik.default = function(m,betahat,betahatsd,v,pilik,FUN="+"){
  sum(log(dens_conv_mixlik(m,betahat,betahatsd,v,pilik,FUN)))
}

#compute the density of the components of the mixture m
#when convoluted with l-components normal mixture with standard deviation s
#or C (scale vector) multiplies scaled (se) student.t l-components mixture with df v
#with mixture proportion pilik
#the density is evaluated at x
#x is an n-vector
#s and pilik are n by l matrices
#v and c are l-vectors
#m is a mixture with k components
#output is a (k*l) by n matrix of densities

#' @title compdens_conv_mixlik
#' @description compute the density of the components of the mixture m
#'     when convoluted with l-components normal mixture with standard
#'     deviation s or C (scale vector) multiplies scaled (se)
#'     student.t l-components mixture with df v with mixture
#'     proportion pilik the density is evaluated at x
#'
#' @param m a mixture with k components
#' @param x an n vector
#' @param s normal mixture of sd(s), n-by-l matrix
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion , n-by-l matrix
#' @param FUN default is "+"
#' @return a (k*l) by n matrix of densities
#' @export
#'
#todo: C is missing for function input
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


#' @title dens_conv_mixlik
#' @description compute density of mixture m convoluted with
#'     l-components normal mixture of sd (s) or student t mixture with
#'     df v with mixture proportion pilik at locations x
#'
#' @param m mixture distribution
#' @param x an n vector
#' @param s normal mixture of sd(s), n-by-l matrix
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion , n-by-l matrix
#' @param FUN default is "+"
#' @export
#'
dens_conv_mixlik = function(m,x,s,v,pilik,FUN="+"){
  UseMethod("dens_conv_mixlik")
}
dens_conv_mixlik.default = function(m,x,s,v,pilik,FUN="+"){
  l=dim(pilik)[2]
  colSums(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik,FUN))
}


#' @title comppostprob_mixlik
#' @description compute the posterior prob that each observation came
#'     from each component of the mixture m,output a k by n vector of
#'     probabilities computed by weighting the component densities by
#'     pi and then normalizing,when likelihood is an l-components
#'     normal mixture or student t mixture with mixture proportion
#'     pilik
#'
#' @param m mixture distribution
#' @param x an n vector
#' @param s normal mixture of sd(s), n-by-l matrix
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion , n-by-l matrix
#' @export
comppostprob_mixlik=function(m,x,s,v,pilik){
  UseMethod("comppostprob_mixlik")
}
#' @export
comppostprob_mixlik.default = function(m,x,s,v,pilik){
  l=dim(pilik)[2]
  k=length(m$pi)
  tmp= (t(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik))/dens_conv_mixlik(m,x,s,v,pilik))
  ismissing = (is.na(x) | apply(is.na(s),1,sum))
  tmp[ismissing,]=rep(m$pi,l)/l
  group=rep(1:k,l)
  return(rowsum(t(tmp),group))
}


#' @title comppostprob_mixlik2
#' @description compute the posterior prob that each observation came
#'     from each component of the mixture m and the likelihood
#'     mixture, output a (k*l) by n vector of probabilities computed
#'     by weighting the component densities by pi and then
#'     normalizing, when likelihood is an l-components normal mixture
#'     or student t mixture with mixture proportion pilik.
#'
#' @param m mixture distribution
#' @param x an n vector
#' @param s normal mixture of sd(s), n-by-l matrix
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion , n-by-l matrix
#' @export
comppostprob_mixlik2=function(m,x,s,v,pilik){
  UseMethod("comppostprob_mixlik2")
}
#' @export
comppostprob_mixlik2.default = function(m,x,s,v,pilik){
  l=dim(pilik)[2]
  k=length(m$pi)
  tmp= (t(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik))/dens_conv_mixlik(m,x,s,v,pilik))
  ismissing = (is.na(x) | apply(is.na(s),1,sum))
  tmp[ismissing,]=rep(m$pi,l)/l
  return(t(tmp))
}



#' @title compcdf_post_mixlik
#' @description evaluate cdf of posterior distribution of beta at c
#' @param m prior on beta, a mixture
#' @param c location of evaluation
#' @param betahat the data
#' @param sebetahat the observed standard errors
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion
#' @return it returns a (k*l) by n matrix
#' @export
compcdf_post_mixlik=function(m,c,betahat,sebetahat,v,pilik){
  UseMethod("compcdf_post_mixlik")
}
#' @export
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


#' @title postmean_mixlik
#' @description output posterior mean for beta for prior mixture m
#' @param m mixture distribution
#' @param betahat the data
#' @param sebetahat the observed standard errors
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion
#' @export
#'
postmean_mixlik = function(m, betahat,sebetahat,v,pilik){
  UseMethod("postmean_mixlik")
}
postmean_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  colSums(comppostprob_mixlik2(m,betahat,sebetahat,v,pilik) * comp_postmean_mixlik(m,betahat,sebetahat,v,pilik))
}


#' @title postmean2_mixlik
#' @description output posterior mean-squared value for beta for prior mixture m,given observations betahat, sebetahat, df v, from l-components mixture likelihood with mixture proportion pilik
#' @param m mixture distribution
#' @param betahat the data
#' @param sebetahat the observed standard errors
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion
#' @export
postmean2_mixlik = function(m, betahat,sebetahat,v,pilik){
  UseMethod("postmean2_mixlik")
}
#' @export
postmean2_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  colSums(comppostprob_mixlik2(m,betahat,sebetahat,v,pilik) * comp_postmean2_mixlik(m,betahat,sebetahat,v,pilik))
}

#' @title postsd_mixlik
#' @description output posterior sd for beta for prior mixture m
#' @param m mixture distribution
#' @param betahat the data
#' @param sebetahat the observed standard errors
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion
#' @export
#'
postsd_mixlik = function(m,betahat,sebetahat,v,pilik){
  UseMethod("postsd_mixlik")
}
postsd_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  sqrt(postmean2_mixlik(m,betahat,sebetahat,v,pilik)-postmean_mixlik(m,betahat,sebetahat,v,pilik)^2)
}


#' @title comp_postmean2_mixlik
#' @description output posterior mean-squared value for beta
#' @param m mixture distribution
#' @param betahat the data
#' @param sebetahat the observed standard errors
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion
#' @export
comp_postmean2_mixlik = function(m,betahat,sebetahat,v,pilik){
  UseMethod("comp_postmean2_mixlik")
}
#' @export
comp_postmean2_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  comp_postsd_mixlik(m,betahat,sebetahat,v,pilik)^2 +
    comp_postmean_mixlik(m,betahat,sebetahat,v,pilik)^2
}


#' @title comp_postmean_mixlik
#' @description output posterior mean for beta for each component of
#'     prior mixture m,given observations betahat, sebetahat, df v
#'     from l-components mixture likelihood with mixture proportion
#'     pilik
#' @param m mixture distribution
#' @param betahat the data
#' @param sebetahat the observed standard errors
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion
#' @export
#'
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


#' @title comp_postsd_mixlik
#'
#' @description output posterior sd for beta for each component of
#'     prior mixture m,given observations betahat, sebetahat, df v
#'     from l-components mixture likelihood with mixture proportion
#'     pilik
#' @param m mixture distribution
#' @param betahat the data
#' @param sebetahat the observed standard errors
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion
#' @export
#'
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

