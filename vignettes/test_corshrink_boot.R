

################  testing the Bootstrap CorShrink models  ########################

dim <- 100
nsamples <- 10
generate_sample <- mvtnorm::rmvnorm(nsamples, rep(0, dim), pop_cov);

sd_boot <- bootcorSE_calc(generate_sample, nboot=50)

cormat <- cor(generate_sample)


generate_sample <- mvtnorm::rmvnorm(nsamples, rep(0, dim), pop_cov);
cov_sample <- cov(generate_sample)
eigen_sample <- eigen(cov_sample, only.values = TRUE)
nsamples <- nsamples

system.time(cov_sample_ML <-  CorShrink::CovShrink(generate_sample, nsamples,
                                                   sd_boot = TRUE,
                                                   type="ML"))

system.time(cov_sample_ML_2 <-  CorShrink::CovShrink(generate_sample, nsamples,
                                                   sd_boot = FALSE,
                                                   type="ML"))

par(mfrow=c(1,2))
image(cov2cor(as.matrix(cov_sample_ML)))
image(cov2cor(as.matrix(cov_sample_ML_2)))

image(cov2cor(cov_sample))
image(cov2cor())
