#' Gene Atlas Human (GenAtHum) Data
#'
#' @description
#' The microarray data were obtained from the GDS DataSet of the Gene
#' Expression Omnibus (GEO) repository in the NCBI archives.
#'
#' @format A list of length 6 containing:
#' \describe{
#'   \item{X}{Design matrix of dimension 158 x 2045.}
#'   \item{y}{Response vector of length 158.}
#'   \item{group}{Group membership code for each predictor.}
#'   \item{gene_code}{Vector of gene codes.}
#'   \item{gene_name}{Vector of gene names.}
#'   \item{groups_name}{Vector of group names.}
#' }
#'
#' @details
#' This data set contains 158 samples, 2045 predictors, and 79 groups,
#' as described in Yuzbasi and Cao (2025).
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/}
#'
#' @examples
#' data(GenAtHum)
#' X <- GenAtHum$X
#' y <- GenAtHum$y
#' group <- GenAtHum$group
#'
#' @docType data
#' @name GenAtHum
NULL