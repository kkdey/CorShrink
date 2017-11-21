
##############################  test  1   ###################################


common_samples_mat <- get(load(file = "../data/common_samples.rda"))
data <- get(load(file = "../data/sample_by_feature_data.rda"))
cormat <- cor(data, use = "pairwise.complete.obs")

