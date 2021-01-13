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
#' @param whichRowTree Which row tree to be replaced? Default is 1 (the first
#'   tree in the \code{rowTree} slot).
#' @param whichColTree Which column tree to be replaced? Default is 1 (the first
#'   tree in the \code{colTree} slot).
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
#' # if rownames are not used in node labels of the tree, provide rowNodeLab
#' use <- changeTree(x = tse, rowTree = treeR2, 
#'                   rowNodeLab = rep(treeR$tip.label, each =2))
#' use
#' 
#' # if rownames are used in node labels of tree, rowNodeLab is not required.
#' 
#' rownames(tse) <- rep(treeR$tip.label, each =2)
#' cse <- changeTree(x = tse, rowTree = treeR2)
#' cse
changeTree <- function(x, 
                       rowTree = NULL, rowNodeLab = NULL,
                       colTree = NULL, colNodeLab = NULL,
                       whichRowTree = 1, whichColTree = 1) {
    
    if (!is.null(rowTree)) {
        out <- .replace_tree(x = x, value = rowTree, 
                             whichTree = whichRowTree, 
                             nodeLab = rowNodeLab, dim = "row")
        if (length(out$drop)) {x <- x[-out$drop, ]}
        x <- BiocGenerics:::replaceSlots(object = x, 
                                         rowTree = out$new_tree,
                                         rowLinks = out$new_links)
    }
    
    if (!is.null(colTree)) {
        out <- .replace_tree(x = x, value = colTree, 
                             whichTree = whichColTree, 
                             nodeLab = colNodeLab, dim = "col")
        if (length(out$drop)) {x <- x[, -out$drop]}
        x <- BiocGenerics:::replaceSlots(object = x, 
                                         colTree = out$new_tree,
                                         colLinks = out$new_links)
    }
    
    return(x)
}

