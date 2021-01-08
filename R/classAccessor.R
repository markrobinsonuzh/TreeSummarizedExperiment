
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
#' \code{\link[SingleCellExperiment:SingleCellExperiment]{SingleCellExperiment}}
#' should work on \strong{TreeSummarizedExperiment}. Additionally, new accessors
#' \code{rowLinks} \code{colLinks}, \code{rowTree} and \code{colTree} accessor
#' function are available for \strong{TreeSummarizedExperiment}.
#' 
#' @param x A TreeSummarizedExperiment object
#' @param i,j The row, column index to subset \code{x}. The arguments of the
#'   subset function \code{[]}
#' @param drop A logical value, TRUE or FALSE. The argument from the subset
#'   function \code{[]}
#' @param value
#' \itemize{
#'   \item{the new rownames or colnames as a \code{character} value. See
#'     \code{\link[BiocGenerics:row_colnames]{BiocGenerics}}.}
#'   \item{A \code{\link[Biostrings:XStringSet-class]{DNAStringSet}}
#'     object or an object coercible to one}
#' }
#' @param whichTree A numeric indicator or name character to specify which tree
#'   in the \code{rowTree} or \code{colTree} to be extracted. The default is to
#'   extract the first tree. If \code{whichTree = NULL}, a list of all trees is
#'   extracted.
#' @param ... The argument from the subset function \code{[]}
#' @param rowNode A vector of nodes that are used to subset rows. One could use
#'   the node number, the node label or the node alias to specify nodes, but not
#'   a mixture of them.
#' @param colNode A vector of nodes that are used to subset columns. One could
#'   use the node number, the node label or the node alias to specify nodes, but
#'   not a mixture of them.
#' @name TreeSummarizedExperiment-accessor
#' @return Elements from \code{TreeSummarizedExperiment}.
#' @seealso \code{\link{TreeSummarizedExperiment}}
#'   \code{\link[SingleCellExperiment:SingleCellExperiment]{SingleCellExperiment}}
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
#' (rowL <- rowLinks(x = toy_tse))
#' # on columns
#' (colL <- colLinks(x = toy_tse))
#'
#'  ## extract the treeData
#' # on rows
#' (rowT <- rowTree(x = toy_tse))
#' # on columns
#' (colT <- colTree(x = toy_tse))
#'
#' # the referenceSeq data
#' refSeq <- DNAStringSetList(one = DNAStringSet(rep("A",nrow(toy_tse))),
#'                            two = DNAStringSet(rep("B",nrow(toy_tse))))
#' referenceSeq(toy_tse) <- refSeq
#' toy_tse
NULL
