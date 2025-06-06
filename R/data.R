#' @title Gene Atlas Human (GenAtHum) Data
#' 
#' @description 
#' The microarray data is obtained from the GDS DataSet of the Gene Expression Omnibus (GEO) repository in the NCBI archives (\url{www.ncbi.nlm.nih.gov/geo})
#' 
#' @format An object of class list of length 6.
#' 
#' @details
#' This data set contains 158 samples with 2045 predictors, and 79 groups as described in Yuzbasi and Cao (2025).
#' 
#' 
#' @param X is a design matrix [158 x 2045]
#' @param y is a response vector [1 X 158]
#' @param group is the code for each group name
#' @param gene_code is the vector of code for each gene
#' @param gene_name is the vector of names for each gene
#' @param groups_name is the vector of names for each groups
#' @examples
#' data(GenAtHum)
#' X <- GenAtHum$X
#' y <- GenAtHum$y
#' group <- GenAtHum$group
"GenAtHum"



