

##################   Subspace comparison analysis  ###################


dim <- 100
nsamples <- 5
pop_cov <- generate_cov(dim);

generate_sample <- mvtnorm::rmvnorm(nsamples, rep(0, dim), pop_cov$Sigma);
cov_sample <- cov(generate_sample)
eigen_sample <- eigen(cov_sample, only.values = TRUE)
nsamples <- nsamples

system.time(cov_sample_ML <-  CorShrink::CovShrink(generate_sample, nsamples, sd_boot = FALSE, type="ML"))
system.time(cov_sample_VEM <-  CorShrink::CovShrink(generate_sample, nsamples, type="VEM"))
system.time(strimmer_sample <- corpcor::cov.shrink(generate_sample))
system.time(glasso_sample_005 <- glasso::glasso(cov_sample, rho = 0.05))
system.time(glasso_sample_05 <- glasso::glasso(cov_sample, rho = 0.5))
system.time(glasso_sample_1 <- glasso::glasso(cov_sample, rho = 1))
system.time(glasso_sample_10 <- glasso::glasso(cov_sample, rho = 10))


library(pracma)
subspace(pop_cov$Sigma, as.matrix(cov_sample_ML))
subspace(pop_cov$Sigma, as.matrix(cov_sample_VEM))
subspace(pop_cov$Sigma, strimmer_sample)
subspace(pop_cov$Sigma, cov_sample)
subspace(pop_cov$Sigma, glasso_sample_005$w)
subspace(pop_cov$Sigma, glasso_sample_05$w)
subspace(pop_cov$Sigma, glasso_sample_1$w)
subspace(pop_cov$Sigma, glasso_sample_10$w)
