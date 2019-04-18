
# ============================================================================
### The documentation of accessor functions
### The corresponding codes are in the file allGenerics.R
# ============================================================================

# -----------------------------------------------------------------------------
### TreeSummarizedExperiment
# -----------------------------------------------------------------------------
#' TreeSummarizedExperiment-accessors
#'
#' All accessor functions that work on
#' \code{\link[SingleCellExperiment]{SingleCellExperiment-class}} should work on
#' \strong{TreeSummarizedExperiment}. Additionally, new accessors \code{rowLink}
#' \code{colLink}, \code{rowTree} and \code{colTree} accessor function are
#' available for \strong{TreeSummarizedExperiment}.
#'
#' @param x A TreeSummarizedExperiment object
#' @param i,j The row, column index to subset \code{x}. The arguments of the
#'   subset function \code{[]}
#' @param drop A logical value, TRUE or FALSE. The argument from the subset
#'   function \code{[]}
#' @param ... The argument from the subset function \code{[]}
#' @name TreeSummarizedExperiment-accessor
#' @return Elements from \code{TreeSummarizedExperiment}.
#' @seealso \code{\link{TreeSummarizedExperiment}}
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment-class}}
#'
#' @author Ruizhu HUANG
#' @examples
#'
#' # the assay table
#' set.seed(1)
#' y <- matrix(rnbinom(300,size=1,mu=10),nrow=10)
#' colnames(y) <- paste(rep(LETTERS[1:3], each = 10), rep(1:10,3), sep = "_")
#' rownames(y) <- tinyTree$tip.label
#'
#' # the row data
#' rowInf <- DataFrame(var1 = sample(letters[1:3], 10, replace = TRUE),
#'                     var2 = sample(c(TRUE, FALSE), 10, replace = TRUE))
#' # the column data
#' colInf <- DataFrame(gg = factor(sample(1:3, 30, replace = TRUE)),
#'                     group = rep(LETTERS[1:3], each = 10))
#'
#' # the tree structure on the rows of assay tables
#' data("tinyTree")
#'
#' # the tree structure on the columns of assay tables
#' sampTree <- ape::rtree(30)
#' sampTree$tip.label <- colnames(y)
#'
#' # create the TreeSummarizedExperiment object
#' rowInf$nodeLab <- rownames(y)
#' colInf$nodeLab <- colnames(y)
#' toy_tse <- TreeSummarizedExperiment(assays = list(y),
#'                                     rowData = rowInf,
#'                                     colData = colInf,
#'                                     rowTree = tinyTree,
#'                                     colTree = sampTree)
#'
#' ## extract the rowData
#' (rowD <- rowData(x = toy_tse))
#'
#' ## extract the colData
#' (colD <- colData(x = toy_tse))
#'
#' ## extract the linkData
#' # on rows
#' (rowL <- rowLink(x = toy_tse))
#' # on columns
#' (colL <- colLink(x = toy_tse))
#'
#'  ## extract the treeData
#' # on rows
#' (rowT <- rowTree(x = toy_tse))
#' # on columns
#' (colT <- colTree(x = toy_tse))
#'
NULL

