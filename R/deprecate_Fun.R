#' @rdname findDescendant
#' @export
findOS <- function(tree,
                   node,
                   only.leaf = TRUE,
                   self.include = FALSE,
                   use.alias = FALSE){
    .Deprecated("findDescendant")
    findDescendant(tree = tree,
                   node = node,
                   only.leaf = only.leaf,
                   self.include = self.include,
                   use.alias = use.alias)
}

#' @rdname convertNode
#' @export
transNode <- function(tree, node, use.alias = FALSE,
                      message = FALSE){
    .Deprecated("convertNode")
    convertNode(tree = tree, node = node, use.alias = use.alias,
                message = message)
}


#' Perform data aggregations based on the available tree structures
#'
#' \code{aggValue} aggregates values on the leaf nodes of a tree to a specific
#' arbitrary level of the tree. The level is specified via the nodes of the
#' tree. Users could decide on which dimension (row or column) and how should
#' the aggregation be performed.
#'
#' @param x A \code{TreeSummarizedExperiment} object.
#' @param rowLevel A numeric (node numbers) or character (node labels) vector.
#'   It provides the level on the tree that data is aggregated to. The
#'   aggregation is on the row dimension. The default is \code{rowLevel = NULL},
#'   and no aggregation is performed.
#' @param rowBlock A column name in the \code{rowData} to separate the
#'   aggregation.
#' @param colLevel A numeric (node numbers) or character (node labels) vector.
#'   It provides the level on the tree that data is aggregated to. The
#'   aggregation is on the column dimension. The default is \code{colLevel =
#'   NULL}, and no aggregation is performed.
#' @param colBlock A column name in the \code{colData} to separate the
#'   aggregation.
#' @param FUN A function to be applied on the aggregation. It's similar to the
#'   \code{FUN} in \code{\link[base]{apply}}.
#' @param assay A integer scalar or string indicating which assay of \code{x} to
#'   use in the aggregation. If NULL, all assay tables are used in aggregation.
#' @param message A logical value. The default is TRUE. If TRUE, it will print
#'   out the running process.
#' @return A \code{TreeSummarizedExperiment} object or a
#'   \code{matrix}. The output has the same class of the input \code{x}.
#' @export
#'
#' @author Ruizhu HUANG
#' @seealso \code{\link{aggTSE}}
aggValue <- function(x, rowLevel = NULL, rowBlock = NULL, 
                     colLevel = NULL, colBlock = NULL, 
                     FUN = sum, assay = NULL,
                     message = FALSE){
    .Deprecated("aggTSE")
    aggTSE(x = x, 
           rowLevel = rowLevel, rowBlock = rowBlock, 
           colLevel = colLevel, colBlock = colBlock, 
           rowFun = FUN, colFun = FUN, 
           rowDataCols = colnames(rowData(x)),
           colDataCols = colnames(colData(x)),
           whichAssay = assay,
           message = message,
           whichRowTree = 1, whichColTree = 1,
           rowFirst = TRUE)
}