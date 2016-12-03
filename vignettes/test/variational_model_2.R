


###############   Model 2  ###########################


estimate_mixprop_2 <- function(betahat, sebetahat, g, prior,
                               null.comp=1,df=NULL,control=list()){
  pi_init=NULL
  k=ncomp(g)
  n = length(betahat)
  control$tol = min(0.1/n,1.e-7) # set convergence criteria to be more stringent for larger samples

  matrix_llik = t(log_compdens_conv(g,betahat,sebetahat,df)) #an n by k matrix
  matrix_llik = matrix_llik - apply(matrix_llik,1, max) #avoid numerical issues by subtracting max of each row
  matrix_lik = exp(matrix_llik)

  fit=do.call("mixVBEM2",args = list(matrix_lik= matrix_lik, prior=prior, pi_init=pi_init, control=control))


  null.loglik = sum(log(matrix_lik[,null.comp]))
  loglik.final =  penloglik(fit$pihat,matrix_lik,1)
  g$pi=pi
  if(controlinput$trace==TRUE){toc()}

  ll <- list(pi= fit$pihat, nu1 = fit$nu1hat,
             nu2 = fit$nu2hat, loglik=loglik.final,
             null.loglik=null.loglik, matrix_lik=matrix_lik)
  return(ll)
}


mixVBEM2 = function(matrix_lik, prior, pi_init = NULL,
                    nu1_init = NULL, nu2_init = NULL,
                    nullweight = 100, control=list()){
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)

  k=ncol(matrix_lik)
  if(is.null(pi_init)){  pi_init = c(nullweight, rep(1, k-1))  }# Use as starting point for pi
  if(is.null(nu1_init)){ nu1_init = rep(1, k)}
  if(is.null(nu2_init)){ nu2_init = rep(1, k)}

  param_init <- c(pi_init, log(nu1_init+1e-15, base=2), log(nu2_init+1e-15, base=2))

  res = squarem(par=param_init,fixptfn=VBfixpoint2, objfn=VBnegpenloglik2,matrix_lik=matrix_lik, prior=prior, control=controlinput)

  pi_new = res$par[1:k]
  nu1_new = res$par[(k+1):(2*k)]
  nu2_new = res$par[(2*k+1):(3*k)]

  return(list(pihat = pi_new/sum(pi_new),
              nu1hat = nu1_new,
              nu2hat = nu2_new,
              B=res$value.objfn,
              niter = res$iter,
              converged=res$convergence,
              post=res$par))
}


VBfixpoint2 = function(param_post, matrix_lik, prior){
  n=nrow(matrix_lik)
  k=ncol(matrix_lik)

  pi_post = param_post[1:k];
  nu1_post = 2^(param_post[(k+1):(2*k)]) - 1e-15
  nu2_post = 2^(param_post[(2*k+1):(3*k)]) - 1e-15

  mat1 <- matrix(exp(rep(digamma(pi_post),n)-rep(digamma(sum(pi_post)),k*n) + 0.5*rep(digamma(nu1_post),n) - 0.5* rep(log(nu2_post),n)),ncol=k,byrow=TRUE)
  mat2 <- exp(- (betahat^2/2)%*% t(nu1_post/nu2_post));

  classprob <- mat1*mat2;
  classprob = classprob/rowSums(classprob) # n by k matrix

  pi_post_new <- prior + colSums(classprob)
  pi_post_new <- pi_post_new/ sum(pi_post_new)

  nu1_new <- nu1_post + 0.5*colSums(mat_normed)

  tmp <- matrix(rep(betahat^2,k), ncol=k, byrow=TRUE)*mat_normed
  nu2_new <- nu2_post + colSums(tmp)

  param_new <- c(pi_post_new, log(nu1_new + 1e-15), log(nu2_new + 1e-15));

  return(param_new)
}

VBnegpenloglik2=function(param_post,matrix_lik,prior){
  return(-VBpenloglik2(param_post,matrix_lik,prior))
}

VBpenloglik2 = function(param_post, matrix_lik, prior){
  n=nrow(matrix_lik)
  k=ncol(matrix_lik)

  pi_post = param_post[1:k];
  nu1_post = exp(param_post[(k+1):(2*k)]) - 1e-15
  nu2_post = exp(param_post[(2*k+1):(3*k)]) - 1e-15

  mat1 <- matrix(exp(rep(digamma(pi_post),n)-rep(digamma(sum(pi_post)),k*n) + 0.5*rep(digamma(nu1_post),n) - 0.5* rep(log(nu2_post),n)),ncol=k,byrow=TRUE)
  mat2 <- exp(- (betahat^2/2)%*% t(nu1_post/nu2_post));

  classprob1 <- mat1*mat2;
  classprob = classprob1/rowSums(classprob1) # n by k matrix

  invgammaloglik = k*(nu1_post + 1)*log(nu2_post) - (nu1_post + 1) *k*digamma(nu1_post) - k * nu1_post

  B= sum(classprob*log(classprob1),na.rm=TRUE) - diriKL(prior,pi_post) - sum(classprob*log(classprob)) + invgammaloglik
  return(B)
}

