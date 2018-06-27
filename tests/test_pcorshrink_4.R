library(CorShrink)
library(huge)
library(corpcor)
DM_toeplitz = function(n,P){
  library("MASS")
  index1=sort(sample(seq(1:n),(n/2)))
  index2=seq(1:n)[-index1]
  Sigmatp=function(P){
    a=array(0,dim=c(P,P))
    for(i in 1:P){
      for(j in 1:P){
        a[i,j]=max(1-0.1*(abs(i-j)),0)
      }
    }
    return(a)
  }
  Sigma = Sigmatp(P)
  data = mvrnorm(n,rep(0,P),Sigma)
  Xtest = data[index2,]
  Xtrain = data[index1,]
  Omega = solve(Sigma)
  return(list(Xtrain = Xtrain, Xtest = Xtest, Sigma = Sigma))
}

n <- 1000
P <- 100
ll <- DM_toeplitz(n=n, P=P)
data <- rbind(ll$Xtrain, ll$Xtest)
Sigma <- ll$Sigma
corSigma <- cov2cor(Sigma)
pcorSigma <- cor2pcor(corSigma)

corrplot(as.matrix(corSigma), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.9, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

corrplot(as.matrix(pcorSigma), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.9, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

pcor2 <- pCorShrinkData(data,
                        reg_type = "lm",
                        glmnet_alpha = 1,
                        ash.control = list(mixcompdist="halfuniform",
                                           control = list(maxiter=1000)))
corrplot(as.matrix(pcor2), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.9, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")


