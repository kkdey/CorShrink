

################   Example simulation Scheme  ##########################

dim <- 100
nsamples <- 50
pop_cov <- generate_cov(dim);
plot(sort(pop_cov$eigen, decreasing = TRUE), type="l")


generate_sample <- mvtnorm::rmvnorm(nsamples, rep(0, dim), pop_cov$Sigma);
cov_sample <- cov(generate_sample)
eigen_sample <- eigen(cov_sample, only.values = TRUE)
nsamples <- nsamples

cov_sample_ML <-  CorShrink::CovShrink(generate_sample, nsamples, type="ML")
cov_sample_VEM <-  CorShrink::CovShrink(generate_sample, nsamples, type="VEM")
cov_sample_VEM2 <-  CorShrink::CovShrink(generate_sample, nsamples, type="VEM2")
strimmer_sample <- corpcor::cov.shrink(generate_sample)
glasso_sample_005 <- glasso::glasso(cov_sample, rho = 0.05)
glasso_sample_05 <- glasso::glasso(cov_sample, rho = 0.5)
glasso_sample_1 <- glasso::glasso(cov_sample, rho = 1)


eigen_sample_ML <- eigen(cov_sample_ML, only.values = TRUE)
eigen_sample_VEM <- eigen(cov_sample_VEM, only.values = TRUE)
eigen_sample_VEM2 <- eigen(cov_sample_VEM2, only.values = TRUE)
eigen_strimmer <- eigen(strimmer_sample, only.values = TRUE)
eigen_glasso_005 <- eigen(glasso_sample_005$w, only.values = TRUE)
eigen_glasso_05 <- eigen(glasso_sample_05$w, only.values = TRUE)
eigen_glasso_1 <- eigen(glasso_sample_1$w, only.values = TRUE)
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
  geom_line(aes(y = sample.cov, colour = "sample.cov"))+
  geom_line(aes(y = pop.cov, colour = "pop.cov"))+
  xlab("Order of eigenvalues (sorted)")+
  ylab("log(Eigenvalues + 1)")+
  scale_colour_manual(values=c("blue", "purple", "magenta", "pink", "orange", "red", "green", "yellow", "black", "grey"))+
  ggtitle(paste0("Eigenvalues distribution n/p=", round(nsamples/dim, 4)," for different shrinkage methods"))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))


#################  image plots for the correlation matrices  ##########################

library(fields)
set.seed(1)
par(mfrow=c(2,2))
cols = gray.colors(100)
image.plot(cov2cor(pop_cov$Sigma), col=cols, nlevel=100)
image.plot(as.matrix(ash_cor_sample), col=cols, nlevel=100)
image.plot(as.matrix(ash_cor_sample2), col=cols, nlevel=100)
image.plot(as.matrix(ash_cor_sample3), col=cols, nlevel=100)
