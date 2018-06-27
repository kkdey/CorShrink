

library(CorShrink)
library(corpcor)
library(glasso)
library(Matrix)
library(psych)
library(corrplot)

n <- 500
P <- 100
NUM_SIM <- 50

diags <- list()
diags[[1]] <- rep(1, 100)
diags[[2]] <- rep(-0.5, 100)
Kinv <- bandSparse(100, k = -(0:1), diag = diags, symm = TRUE)
K <- solve(Kinv)
corSigma <- cov2cor(K)
data <- MASS::mvrnorm(n,rep(0,P),corSigma)
pcorSigma <- cor2pcor(corSigma)

corrplot(as.matrix(pcorSigma), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.9, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

corrplot(as.matrix(corSigma), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.9, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")



pcor2 <- pCorShrinkData(data,
                        reg_type = "glmnet",
                        glmnet_alpha = 1,
                        ash.control = list(mixcompdist="normal", control = list(maxiter=1000)))
corrplot(as.matrix(pcor2), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.9, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

