

##############   Variational EM on CorShrink   #######################


###  We perform adaptive shrinkage using Variational EM on the correlation
#########  matrix.


mat <- generate_cov(dim=100);
Sigma <- mat$Sigma;
corSigma <- cov2cor(Sigma);

nsamples = 50;
cor_table <- reshape2::melt(corSigma);
cor_table_non_diag <- cor_table[which(cor_table[,1] !=cor_table[,2]),];

cor <- ash_cor (corSigma, nsamples = 500, tol=1e-06)


###################  correlation table   ###########################


cor_table_non_diag.val <- cor_table_non_diag[,3];
cor_table_non_diag.val[which(cor_table_non_diag.val==1)]=(1- 1e-7);
#cor_table_non_diag.val[which(cor_table_non_diag.val==0)]=(1e-7);


#########  extracting the Fisher Z-scores and the variances   ##############


cor_transform_mean_vec=0.5*log((1+cor_table_non_diag.val)/(1-cor_table_non_diag.val))
cor_transform_sd_vec=rep(sqrt(1/(nsamples-3)), dim(cor_table_non_diag)[1]);
options(warn=-1)


#########  setting up the sigma values   #######################

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

betahat <- cor_transform_mean_vec
sebetahat <- cor_transform_sd_vec
completeobs <- which(!is.na(betahat));
gridmult = sqrt(2);
mixsd = autoselect.mixsd(betahat[completeobs]-mean(betahat[completeobs]),sebetahat[completeobs],gridmult)
mixsd = c(0,mixsd)
null.comp = which.min(mixsd) #which component is the "null"

k = length(mixsd)
prior = "nullbiased"
nullweight=10
null.comp=1
prior = setprior(prior,k,nullweight,null.comp)
randomstart = FALSE
n = length(betahat)
normalize <- function(x) { return(x/sum(x))}
pi = initpi(k,n,null.comp,randomstart)


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

image.plot(as.matrix(ash_cor_PD2), col=two.colors(n=256, start="darkgreen", end="red", middle="white",
                                                 alpha=1.0))
image.plot(as.matrix(corSigma), col=two.colors(n=256, start="darkgreen", end="red", middle="white",
                                               alpha=1.0))
