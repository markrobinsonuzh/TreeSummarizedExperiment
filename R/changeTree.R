#' Change the row or column tree 
#' 
#' \code{changeTree} changes a row or column tree in a
#' \code{TreeSummarizedExperiment} object.
#' 
#' @param x A TreeSummarizedExperiment object
#' @param rowTree A phylo object. A new row tree.
#' @param rowNodeLab A character string. It provides the labels of nodes that
#'   the rows of assays tables corresponding to. If NULL (default), the row
#'   names of the assays tables are used.
#' @param colTree A phylo object. A new column tree.
#' @param colNodeLab A character string. It provides the labels of nodes that
#'   the columns of assays tables corresponding to. If NULL (default), the
#'   column names of the assays tables are used.
#' @return A TreeSummarizedExperiment object
#' @importFrom stats setNames
#' @export
#' @author Ruizhu Huang
#' @examples 
#' 
#' library(ape)
#' set.seed(1)
#' treeR <- ape::rtree(10)
#' 
#' # the count table
#' count <- matrix(rpois(160, 50), nrow = 20)
#' rownames(count) <- paste0("entity", 1:20)
#' colnames(count) <- paste("sample", 1:8, sep = "_")
#' # The sample information
#' sampC <- data.frame(condition = rep(c("control", "trt"), 
#'                                     each = 4),
#'                     gender = sample(x = 1:2, size = 8, 
#'                                     replace = TRUE))
#' rownames(sampC) <- colnames(count)
#' # build a TreeSummarizedExperiment object
#' tse <- TreeSummarizedExperiment(assays = list(count),
#'                                 colData = sampC,
#'                                 rowTree = treeR,
#'                                 rowNodeLab = rep(treeR$tip.label, each =2))
#' 
#' treeR2 <- drop.tip(phy = treeR, tip = c("t10", "t9", "t8"))
#' 
#' use <- changeTree(x = tse, rowTree = treeR2)
#' use

changeTree <- function(x, 
                       rowTree = NULL, rowNodeLab = NULL,
                       colTree = NULL, colNodeLab = NULL) {
    
    if (!is.null(rowTree)) {
        rlab <- c(rowTree$tip.label, rowTree$node.label)
        if (is.null(rowNodeLab)) {
            rowNodeLab <- intersect(rownames(x), rlab)
            missLab <- setdiff(rownames(x), rlab)
            nm <- length(missLab)
            if (nm) {
                warning(nm, " rows are removed due to mismatch")}
            x <- x[rowNodeLab, ]
        }
        nrLink <- .updateLinks(newTree = rowTree, newLab = rowNodeLab)
        rownames(nrLink) <- rownames(x)
        x <- BiocGenerics:::replaceSlots(x,
                                         rowLinks = nrLink,
                                         rowTree = list(phylo = rowTree))
    }
    
    if (!is.null(colTree)) {
        clab <- c(colTree$tip.label, colTree$node.label)
        if (is.null(colNodeLab)) {
            colNodeLab <- intersect(colnames(x), clab)
            missLab <- setdiff(colnames(x), clab)
            nm <- length(missLab)
            if (nm) {
                warning(nm, " columns are removed due to mismatch")}
            x <- x[, colNodeLab]
        }
        ncLink <- .updateLinks(newTree = colTree, newLab = colNodeLab)
        rownames(ncLink) <- colnames(x)
        x <- BiocGenerics:::replaceSlots(x,
                                         colLinks = ncLink,
                                         colTree = list(phylo = colTree))
    }
    
    
    return(x)
}

.updateLinks <- function(newTree, newLab) {
    num <- convertNode(tree = newTree, node = newLab)
    lab <- convertNode(tree = newTree, node = num)
    alias <- convertNode(tree = newTree,
                       node = num, 
                       use.alias = TRUE)
    ind <- isLeaf(tree = newTree, node = num)
    LinkDataFrame(nodeLab = lab, nodeLab_alias = alias, 
                  nodeNum = num, isLeaf = ind)
}
