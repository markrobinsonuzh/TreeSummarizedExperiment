#' Combine TSEs by rows or columns
#'
#'  
#' \code{rbind} and \code{cbind} take one or more
#' \code{TreeSummarizedExperiment} objects and combine them by columns or rows,
#' respectively.
#' 
#' @param ... One or more \code{TreeSummarizedExperiment} objects.
#' @param deparse.level See \code{\link[base]{cbind}}
#' @name TreeSummarizedExperiment-combine
#' @return A TreeSummarizedExperiment object
#' 
#' @author Ruizhu Huang
#' @examples 
#' 
#' # rbind works : 
#' # a) TSE without rowTree and without colTree
#' # b) TSE with rowTree but without colTree
#' # c) TSE without rowTree but with colTree
#' # d) TSE with rowTree & colTree
#' 
#' set.seed(1)
#' # a) 
#' (tse_a <- makeTSE(include.colTree = FALSE))
#' (tse_b <- makeTSE(include.colTree = FALSE))
#' 
#' # b) 
#' (tse_c <- makeTSE(include.rowTree = FALSE))
#' (tse_d <- makeTSE(include.rowTree = FALSE))
#' 
#' rbind(tse_a, tse_b)
#' cbind(tse_c, tse_d)

NULL