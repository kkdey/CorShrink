
###########  Test pCorShrink  ##############################

require(graphics)
library(igraph)
library(MASS)
library(Matrix)

make.data <- function(Sigma.half, n, p, seed){
  set.seed(seed)

  X = matrix(rnorm(n*p),n,p)%*%Sigma.half
  return(X)
}

band.mat <- function(a, p, K=1, permu=c(1:p)){
  ones = rep(1,p)
  Omega0 = a*ones%*%t(ones)
  diag(Omega0) = rep(1,p)
  Omega = 1*band(Omega0,-K,K)
  Sigma = qr.solve(Omega)
  Sigma = Sigma*(abs(Sigma)>1e-4)
  Sigma.half=chol(Sigma)
  Sigma.half = Sigma.half*(abs(Sigma.half)>1e-4)
  Sigma = Sigma[permu,permu]
  Omega = Omega[permu,permu]
  Sigma.half = Sigma.half[permu,permu]
  obj = list(Sigma=Sigma, Omega = Omega, Sigma.half = Sigma.half)
}

# simulation parameters
n = 200
p = 50
obj = band.mat(a=0.5, p, K = 1)

library(corrplot)

pcor1 <- - cov2cor(obj$Omega)
diag(pcor1) <- 1
col2 <- c("blue", "white", "red")
corrplot(as.matrix(pcor1), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.9, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

Sig.half = obj$Sigma.half
Ome.true = obj$Omega
X.mat = make.data(Sig.half, n, p, seed = 1000)

pcor2 <- pCorShrinkData(X.mat, ash.control = list(mixcompdist="normal", control = list(maxiter=1000)))
corrplot(as.matrix(pcor2), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.9, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

