

library(MASS)
library(Matrix)
library(glasso)
library(CorShrink)
library(corpcor)

P <- 100
block <- 10
mat <- 0.3*diag(1,block) + 0.7*rep(1,block) %*% t(rep(1, block))
Sigma <-   bdiag(mat, mat, mat, mat)
set.seed(100)
P <- dim(Sigma)[1]
n=100
data <- MASS::mvrnorm(n,rep(0,P),Sigma)

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


