########   estimate CorShrink  model  2   ###############################

mat <- generate_cov(dim=100);
Sigma <- mat$Sigma;
corSigma <- cov2cor(Sigma);

nsamples = 50;
cor_table <- reshape2::melt(corSigma);
cor_table_non_diag <- cor_table[which(cor_table[,1] !=cor_table[,2]),];


###################  correlation table   ###########################


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
prior = "nullbiased"
nullweight=10
null.comp=1
prior = setprior(prior,k,nullweight,null.comp)
randomstart = FALSE
n = length(betahat)
normalize <- function(x) { return(x/sum(x))}
pi = initpi(k,n,null.comp,randomstart)


pi_init = NULL
nu1_init = NULL
nu2_init = NULL
g=normalmix(pi,rep(0,k),mixsd)

out <- estimate_mixprop_2(betahat, sebetahat, g, prior)

ghat = g;
ghat$pi = out$pi
ghat$mean = rep(0, k)
ghat$sd = sqrt(out$nu2/(out$nu1 - 1))

postMean.fit <- postmean(ghat,betahat[completeobs],sebetahat[completeobs], v=NULL)

ash_cor_vec=(exp(2*postMean.fit)-1)/(exp(2*postMean.fit)+1);

newdata.table <- cor_table_non_diag;
newdata.table[,3] <- ash_cor_vec;
ash_cor_only <- reshape2::dcast(newdata.table, Var1~Var2, value.var = "value")[,-1];
ash_cor_only[is.na(ash_cor_only)]=1;
pd_completion <- Matrix::nearPD(as.matrix(ash_cor_only), conv.tol=1e-06);
ash_cor_PD <- sweep(pd_completion$mat,diag(as.matrix(pd_completion$mat)), MARGIN=1,"/")

image.plot(as.matrix(ash_cor_PD), col=two.colors(n=256, start="darkgreen", end="red", middle="white",
                                                 alpha=1.0))
image.plot(as.matrix(corSigma), col=two.colors(n=256, start="darkgreen", end="red", middle="white",
                                               alpha=1.0))


