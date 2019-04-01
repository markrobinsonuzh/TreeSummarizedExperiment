
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
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} should work on
#' \strong{TreeSummarizedExperiment}. Additionally, two new \code{linkData} and
#' \code{treeData} accessor function are available for
#' \strong{TreeSummarizedExperiment}.
#'
#' @param x A TreeSummarizedExperiment object
#' @param onRow A logical(1) value, TRUE or FALSE. The default is TRUE. It
#'   indicates the returned \code{treeData} or \code{linkData} is on rows (TRUE)
#'   or columns (FALSE) of \code{assays} tables.
#'
#' @name TreeSummarizedExperiment-accessor
#' @return Elements from \code{TreeSummarizedExperiment}.
#' @seealso \code{\link{TreeSummarizedExperiment}}
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
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
#' (rowL <- linkData(x = toy_tse, onRow = TRUE))
#' # on columns
#' (colL <- linkData(x = toy_tse, onRow = FALSE))
#'
#'  ## extract the treeData
#' # on rows
#' (rowT <- treeData(x = toy_tse, onRow = TRUE))
#' # on columns
#' (colT <- treeData(x = toy_tse, onRow = FALSE))
#'
NULL

# -----------------------------------------------------------------------------
### LinkDataFrame
# -----------------------------------------------------------------------------
#' LinkDataFrame-accessor
#'
#' This page lists the accesssor, coercion, subsetting and combining functions
#' for the class \strong{LinkDataFrame}.
#'
#' @param x A \code{LinkDataFrame} object
#' @param i,j The row, column index to subset \code{x}.
#' @param name The name of the column.
#' @param value The value to replace.
#' @name LinkDataFrame-accessor
#' @return Depends on the functions
#'
#' @author Ruizhu HUANG
#' @seealso \code{\link{LinkDataFrame-constructor}}
#'   \code{\link{LinkDataFrame-class}} \code{\link[S4Vectors]{DataFrame-class}}
#' @examples
#'
#' left <- DataFrame(left1 = 1:5, left2 = letters[1:5])
#' right <- DataFrame(right1 = sample(letters[1:3], 5, replace = TRUE),
#'                   right2 = sample(c(TRUE, FALSE), 5, replace = TRUE),
#'                   right3 = 11:15)
#'
#' ld <- LinkDataFrame(LinkData = left, right)
#'
#' ## subset
#' # by rows
#' (ld1 <- ld[1:3, ])
#' # by columns (only on the right side of the vertical line)
#' (ld2 <- ld[, c(1, 3)])
#' (ld3 <- ld[, c("right1", "right2")])
#'
#' ## add a new column (only on the right side of the vertical line)
#' ld4 <- ld
#' ld4$right4 <- 11:15
#' ld4
#'
#' ## coercion
#' (ld5 <- as.data.frame(ld))
#'
#'
NULL
