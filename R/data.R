#' Samples by features data matrix with missing values
#' A matrix with 544 samples (persons) and 53 features (tissues) of GTEx expression data
#' for gene ENSG00000166819 (https://www.gtexportal.org/home/).
#' @format A matrix with 544 rows and 53 columns, containing missing values coded as NA
#'
"sample_by_feature_data"


#' An example matrix of pairwise correlations
#' A matrix of pairwise correlations computed from the sample_by_feature_data.
#' @format A  symmetric matrix of pairwise correlations with 53 rows and 53 columns.
#' Each entry is a pairwise correlations from the 53 features (tissues) in
#' sample_by_feature_data.
#'
"pairwise_corr_matrix"


#' A matrix of common samples between tissues in sample_by_feature_data
#' Calculates the total non-NA expressions for a pair of tissues (of 53 tissues)
#' and records that in the cell corresponding to each tissue pair
#' @format A matrix of integers specifying number of common persons contributing a
#' pair of tissues for 53 tissues. The matrix is symmetric with 53 rows and 53 columns.
#'
"common_samples"
