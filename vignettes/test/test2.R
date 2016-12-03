

################   Example simulation Scheme  ##########################

dim <- 100
nsamples <- 50
pop_cov <- generate_cov(dim);
plot(sort(pop_cov$eigen, decreasing = TRUE), type="l")


generate_sample <- mvtnorm::rmvnorm(nsamples, rep(0, dim), pop_cov$Sigma);
cov_sample <- cov(generate_sample)
eigen_sample <- eigen(cov_sample, only.values = TRUE)
plot(sort(as.numeric(eigen_sample$values), decreasing = TRUE), type="l")
cor_sample <- cov2cor(cov_sample)


cormat <- cor_sample
nsamples <- nsamples

ash_cor_sample <- ash_cor(cor_sample, nsamples)$ash_cor_PD;
ash_cov_sample <- diag(sqrt(diag(cov_sample)))%*%ash_cor_sample%*%diag(sqrt(diag(cov_sample)))
eigen_ash_sample <- eigen(ash_cov_sample, only.values = TRUE)

ash_cor_sample2 <- ash_cor2(cor_sample, nsamples)$ash_cor_PD;
ash_cov_sample2 <- diag(sqrt(diag(cov_sample)))%*%ash_cor_sample2%*%diag(sqrt(diag(cov_sample)))
eigen_ash_sample2 <- eigen(ash_cov_sample2, only.values = TRUE)

ash_cor_sample3 <- ash_cor3(cor_sample, nsamples)$ash_cor_PD;
ash_cov_sample3 <- diag(sqrt(diag(cov_sample)))%*%ash_cor_sample3%*%diag(sqrt(diag(cov_sample)))
eigen_ash_sample3 <- eigen(ash_cov_sample3, only.values = TRUE)

ash_cov_master_sample <- ash_cov_master(generate_sample);
shrinkage <- ash_cov_master_sample$shrink_intensity
ash_cov_master <- ash_cov_master_sample$ash_cov_ledoit_wolf;
eigen_ash_cov_master <- eigen(ash_cov_master);


library(ggplot2)

eigendata <- data.frame(
  eigen_order = 1:dim,
  ash_cov2 = sort(log(as.numeric(eigen_ash_sample2$values)+1),  decreasing=TRUE),
  ash_cov3 = sort(log(as.numeric(eigen_ash_sample3$values)+1),  decreasing=TRUE),
  ash_cov = sort(log(as.numeric(eigen_ash_sample$values)+1),
                 decreasing = TRUE),
  sample_cov = sort(log(as.numeric(eigen_sample$values)+1),
                    decreasing=TRUE),
#  shafer_strimmer_cov= sort(log(as.numeric(shafer_eigen$values)+1),                   decreasing = TRUE),
  pop_cov = sort(log(as.numeric(pop_cov$eigen)+1), decreasing=TRUE))

colnames(eigendata) <- c("eigenorder",
                         "ash_cov2",
                         "ash_cov3",
                         "ash_cov",
                         "sample_cov",
                         "pop_cov")

library(ggplot2)
ggplot(eigendata, aes(eigenorder)) +
  geom_line(aes(y = ash_cov2, colour = "ash cov2")) +
  geom_line(aes(y = ash_cov3, colour = "ash cov3")) +
  geom_line(aes(y = ash_cov, colour = "ash cov"))+
  geom_line(aes(y = sample_cov, colour = "sample cov"))+
  geom_line(aes(y = pop_cov, colour = "pop cov"))+
  xlab("Order of eigenvalues (sorted)")+
  ylab("log(Eigenvalues + 1)")+
  scale_colour_manual(values=c("blue", "grey", "red", "green", "orange"))+
  ggtitle(paste0("Eigenvalues distribution n/p=", round(nsamples/dim, 4),"\n for different shrinkage methods"))+
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
