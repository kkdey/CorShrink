

################   Example simulation Scheme  ##########################

dim <- 100
nsamples <- 5
pop_cov <- generate_cov(dim);
plot(sort(pop_cov$eigen, decreasing = TRUE), type="l")


generate_sample <- mvtnorm::rmvnorm(nsamples, rep(0, dim), pop_cov$Sigma);
cov_sample <- cov(generate_sample)
eigen_sample <- eigen(cov_sample, only.values = TRUE)
nsamples <- nsamples

system.time(cov_sample_ML <-  CorShrink::CovShrink(generate_sample, nsamples, type="ML"))
system.time(cov_sample_VEM <-  CorShrink::CovShrink(generate_sample, nsamples, type="VEM"))
system.time(cov_sample_VEM2 <-  CorShrink::CovShrink(generate_sample, nsamples, type="VEM2"))
system.time(strimmer_sample <- corpcor::cov.shrink(generate_sample))
system.time(glasso_sample_005 <- glasso::glasso(cov_sample, rho = 0.05))
system.time(glasso_sample_05 <- glasso::glasso(cov_sample, rho = 0.5))
system.time(glasso_sample_1 <- glasso::glasso(cov_sample, rho = 1))
system.time(glasso_sample_10 <- glasso::glasso(cov_sample, rho = 10))



eigen_sample_ML <- eigen(cov_sample_ML, only.values = TRUE)
eigeneigen_sample_VEM <- eigen(cov_sample_VEM, only.values = TRUE)
eigen_sample_VEM2 <- eigen(cov_sample_VEM2, only.values = TRUE)
eigen_strimmer <- eigen(strimmer_sample, only.values = TRUE)
eigen_glasso_005 <- eigen(glasso_sample_005$w, only.values = TRUE)
eigen_glasso_05 <- eigen(glasso_sample_05$w, only.values = TRUE)
eigen_glasso_1 <- eigen(glasso_sample_1$w, only.values = TRUE)
eigen_glasso_10 <- eigen(glasso_sample_10$w, only.values = TRUE)
eigen_sample <- eigen(cov_sample, only.values = TRUE)


library(ggplot2)

eigendata <- data.frame(
  eigen_order = 1:dim,
  covshrink.ML = sort(log(as.numeric(eigen_sample_ML$values)+1),  decreasing=TRUE),
  covshrink.VEM = sort(log(as.numeric(eigen_sample_VEM$values)+1),  decreasing=TRUE),
  covshrink.VEM2 = sort(log(as.numeric(eigen_sample_VEM2$values)+1),
                 decreasing = TRUE),
  covshrink.strimmer = sort(log(as.numeric(eigen_strimmer$values)+1),  decreasing=TRUE),
  glasso.cov.005 = sort(log(as.numeric(eigen_glasso_005$values)+1),  decreasing=TRUE),
  glasso.cov.05 = sort(log(as.numeric(eigen_glasso_05$values)+1),  decreasing=TRUE),
  glasso.cov.1 = sort(log(as.numeric(eigen_glasso_1$values)+1),  decreasing=TRUE),
  glasso.cov.10 = sort(log(as.numeric(eigen_glasso_10$values)+1),  decreasing=TRUE),
  sample_cov = sort(log(as.numeric(eigen_sample$values)+1),
                    decreasing=TRUE),
  pop_cov = sort(log(as.numeric(pop_cov$eigen)+1), decreasing=TRUE))

colnames(eigendata) <- c("eigenorder",
                         "covshrinkML",
                         "covshrinkVEM",
                         "covshrinkVEM2",
                         "cov.strimmer",
                         "cov.glasso.005",
                         "cov.glasso.05",
                         "cov.glasso.1",
                         "cov.glasso.10",
                         "sample.cov",
                         "pop.cov")

library(ggplot2)
ggplot(eigendata, aes(eigenorder)) +
  geom_line(aes(y = covshrinkML, colour = "covshrinkML")) +
  geom_line(aes(y = covshrinkVEM, colour = "covshrinkVEM")) +
  geom_line(aes(y = covshrinkVEM2, colour = "covshrinkVEM2"))+
  geom_line(aes(y = cov.strimmer, colour = "cov.strimmer"))+
  geom_line(aes(y = cov.glasso.005, colour = "cov.glasso.rho.0.05"))+
  geom_line(aes(y = cov.glasso.05, colour = "cov.glasso.rho.0.5"))+
  geom_line(aes(y = cov.glasso.1, colour = "cov.glasso.rho.1"))+
  geom_line(aes(y = cov.glasso.10, colour = "cov.glasso.rho.10"))+
  geom_line(aes(y = sample.cov, colour = "sample.cov"))+
  geom_line(aes(y = pop.cov, colour = "pop.cov"))+
  xlab("Order of eigenvalues (sorted)")+
  ylab("log(Eigenvalues + 1)")+
  scale_colour_manual(values=c("blue", "violet", "purple", "magenta", "pink", "orange", "red", "green", "black", "grey"))+
  ggtitle(paste0("n/p=", round(nsamples/dim, 4),""))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))


eigenval_sample_ML <- eigen(cov_sample_ML)$values
eigenval_sample_VEM <- eigen(cov_sample_VEM)$values
eigenval_sample_VEM2 <- eigen(cov_sample_VEM2)$values
eigenval_strimmer <- eigen(strimmer_sample)$values
eigenval_glasso_005 <- eigen(glasso_sample_005$w)$values
eigenval_glasso_05 <- eigen(glasso_sample_05$w)$values
eigenval_glasso_1 <- eigen(glasso_sample_1$w)$values
eigenval_sample <- eigen(cov_sample)$values
popval <- eigen(pop_cov$Sigma)$values

eigenvec_sample_ML <- eigen(cov_sample_ML)$vectors[,order(eigenval_sample_ML, decreasing=TRUE)]
eigenvec_sample_VEM <- eigen(cov_sample_VEM)$vectors[,order(eigenval_sample_VEM, decreasing=TRUE)]
eigenvec_sample_VEM2 <- eigen(cov_sample_VEM2)$vectors[,order(eigenval_sample_VEM2, decreasing=TRUE)]
eigenvec_strimmer <- eigen(strimmer_sample)$vectors[,order(eigenval_strimmer, decreasing=TRUE)]
eigenvec_glasso_005 <- eigen(glasso_sample_005$w)$vectors[,order(eigenval_glasso_005, decreasing=TRUE)]
eigenvec_glasso_05 <- eigen(glasso_sample_05$w)$vectors[,order(eigenval_glasso_05, decreasing=TRUE)]
eigenvec_glasso_1 <- eigen(glasso_sample_1$w)$vectors[,order(eigenval_glasso_1, decreasing=TRUE)]
eigenvec_sample <- eigen(cov_sample)$vectors[,order(eigenval_sample, decreasing=TRUE)]
popvec <- eigen(pop_cov$Sigma)$vectors[,order(popval, decreasing=TRUE)]

sum(abs(popvec[,1:5] - eigenvec_sample_ML[,1:5]))/(5*dim)
sum(abs(popvec[,1:5] - eigenvec_glasso_005[,1:5]))/(5*dim)
sum(abs(popvec[,1:5] - eigenvec_glasso_05[,1:5]))/(5*dim)
sum(abs(popvec[,1:5] - eigenvec_glasso_1[,1:5]))/(5*dim)
sum(abs(popvec[,1:5] - eigenvec_strimmer[,1:5]))/(5*dim)
sum(abs(popvec[,1:5] - eigenvec_sample[,1:5]))/(5*dim)
sum(abs(popvec[,1:5] - eigenvec_sample_VEM[,1:5]))/(5*dim)
sum(abs(popvec[,1:5] - eigenvec_sample_VEM2[,1:5]))/(5*dim)

#################  image plots for the correlation matrices  ##########################

library(fields)
set.seed(1)
par(mfrow=c(2,2))
cols = rainbow(100)
image.plot(cov2cor(pop_cov$Sigma), col=cols, nlevel=100)
image.plot(cov2cor(as.matrix(cov_sample_ML)), col=cols, nlevel=100)
image.plot(cov2cor(as.matrix(cov_sample_VEM)), col=cols, nlevel=100)
image.plot(cov2cor(as.matrix(cov_sample_VEM2)), col=cols, nlevel=100)

install.packages("cape")
myImagePlot(cov2cor(pop_cov$Sigma))
myImagePlot(cov2cor(as.matrix(cov_sample_ML)))
myImagePlot(as.matrix(strimmer_sample))
myImagePlot(as.matrix(cov_sample_ML))
myImagePlot(cov2cor(as.matrix(cov_sample_VEM2)))
myImagePlot(cov2cor(as.matrix(glasso_sample_005$w)))

pop_cor <- cov2cor(pop_cov$Sigma)
qq <- which(pop_cor == min(pop_cor), arr.ind=T)
cor_sample_VEM2 <- cov2cor(strimmer_sample)
cor_sample_VEM2[qq[1], qq[2]] <- min(pop_cor)
cor_sample_VEM2[qq[2], qq[1]] <- min(pop_cor)
myImagePlot(pop_cor)
myImagePlot(cor_sample_VEM2)

par(mfrow=c(3,2))
cor_sample_ML_vec <- as.vector(cov2cor(as.matrix(cov_sample_ML)))
cor_sample_ML_vec <- cor_sample_ML_vec[cor_sample_ML_vec < 1]
#plot(cor_sample_ML_vec, type="l", col="red")
plot(density(cor_sample_ML_vec), type="l", col="red", xlim=c(-1,1), xlab="correlation", main="CorShrink-ML")

cor_sample_VEM_vec <- as.vector(cov2cor(as.matrix(cov_sample_VEM)))
cor_sample_VEM_vec <- cor_sample_VEM_vec[cor_sample_VEM_vec < 1]
#plot(cor_sample_ML_vec, type="l", col="red")
plot(density(cor_sample_VEM_vec), type="l", col="red", xlim=c(-1,1), xlab="correlation", main="CorShrink-VEM")

cor_sample_VEM2_vec <- as.vector(cov2cor(as.matrix(cov_sample_VEM2)))
cor_sample_VEM2_vec <- cor_sample_VEM2_vec[cor_sample_VEM2_vec < 1]
#plot(cor_sample_ML_vec, type="l", col="red")
plot(density(cor_sample_VEM2_vec), type="l", col="red", xlim=c(-1,1), xlab="correlation", main="CorShrink-VEM2")

strimmer_cor_vec <- as.vector(corpcor::cor.shrink(generate_sample))
strimmer_cor_vec <- strimmer_cor_vec[strimmer_cor_vec < 1]
#plot(strimmer_cor_vec, type="l", col="green")
plot(density(strimmer_cor_vec), type="l", col="green", xlim=c(-1,1), xlab="correlation", main="Shafer=Strimmer")

glasso_cor_vec <- as.vector(cov2cor(as.matrix(glasso_sample_1$w)))
glasso_cor_vec <- glasso_cor_vec[glasso_cor_vec < 1]
#plot(strimmer_cor_vec, type="l", col="green")
plot(density(glasso_cor_vec), type="l", col="green", xlim=c(-1,1), xlab="correlation", main="GLASSO (rho=1)")

pop_cor_vec <- as.vector(cov2cor(pop_cov$Sigma))
pop_cor_vec  <- pop_cor_vec [pop_cor_vec  < 1]
#plot(pop_cor_vec , type="l", col="blue")
plot(density(pop_cor_vec), type="l", col="blue", xlim=c(-1,1), xlab="correlation", main="Pop Corr")





