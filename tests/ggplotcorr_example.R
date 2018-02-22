

###########   ggplotcorr   #####################

library(ggcorrplot)
data("pairwise_corr_matrix")

ggcorrplot(pairwise_corr_matrix)

data("pairwise_corr_matrix")
data("common_samples")
out <- CorShrinkMatrix(pairwise_corr_matrix, common_samples, image_corshrink  = TRUE)

ggcorrplot(out$ash_cor_PD)
ggcorrplot(pairwise_corr_matrix)




data("pairwise_corr_matrix")
data("common_samples")
out <- CorShrinkMatrix(pairwise_corr_matrix, common_samples,
                       image_original = FALSE, image_corshrink  = TRUE)
