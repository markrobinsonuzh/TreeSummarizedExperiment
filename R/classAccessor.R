
# -----------------------------------------------------------------------------
### Documentation of the accessor function
### Function code is in the file allGenerics.R
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
### treeSummarizedExperiment
# -----------------------------------------------------------------------------
#' Accessor functions for treeSummarizedExperiment
#'
#' Accessor functions to extract different elements from
#' \strong{treeSummarizedExperiment} object.
#'
#' @param x A treeSummarizedExperiment object
#' @param ... For assay, ... contains \code{use.nodeLab}, which is forwarded to
#'   assays. For rowData, arguments passed through ... are forwarded to mcols.
#' @param use.nodeLab A logical(1), indicating whether the rownames of assay
#'   elements should use node labels (the column \code{nodeLab} in
#'   \code{linkData} if there is not duplicated values; otherwise the column
#'   \code{nodeLab_alias} in \code{linkData} is used.)
#' @param withDimnames A logical(1), indicating whether dimnames should be
#'   applied to extracted assay elements. Setting withDimnames=FALSE increases
#'   the speed and memory efficiency with which assays are extracted.
#'   withDimnames=TRUE in the getter assays<- allows efficient complex
#'   assignments (e.g., updating names of assays, names(assays(x,
#'   withDimnames=FALSE)) = ... is more efficient than names(assays(x)) = ...);
#'   it does not influence actual assignment of dimnames to assays.
#' @param value An object of a class specified in the S4 method signature or as
#'   outlined in 'Details'.
#' @param i,j The subscripts that can act to subset the rows and columns of
#'   \code{x}, that is the matrix elements of assays.
#' @name treeSummarizedExperiment-accessor
#' @return Elements from \code{treeSummarizedExperiment}.
#' @author Ruizhu HUANG
#' @seealso \code{\link{treeSummarizedExperiment}}
#'   \code{\link{treeSummarizedExperiment-accessor}}
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @examples
#' library(S4Vectors)
#' set.seed(1)
#' y <- matrix(rnbinom(300,size=1,mu=10),nrow=10)
#' colnames(y) <- paste(rep(LETTERS[1:3], each = 10), rep(1:10,3), sep = "_")
#' rownames(y) <- tinyTree$tip.label
#'
#' rowInf <- DataFrame(nodeLab = rownames(y),
#'                     var1 = sample(letters[1:3], 10, replace = TRUE),
#'                     var2 = sample(c(TRUE, FALSE), 10, replace = TRUE))
#' colInf <- DataFrame(gg = factor(sample(1:3, 30, replace = TRUE)),
#'                     group = rep(LETTERS[1:3], each = 10))
#' toy_tse <- treeSummarizedExperiment(tree = tinyTree, rowData = rowInf,
#'                                     colData = colInf,
#'                                     assays = list(y, (2*y), 3*y))
#' rowData(toy_tse)
#' colData(toy_tse)
#' metadata(toy_tse)
#' linkData(toy_tse)
#' assays(toy_tse)
NULL


# @param internal TRUE or FALSE. Only for \code{rowData}. If TRUE, the columns
# with \code{int_rowData} class are visible; otherwise, they would be hiden.
# These columns are usually result tables that users obtain from their
# customized analysis and  are written back to the
# \code{treeSummarizedExperiment} object that stores the original data.
