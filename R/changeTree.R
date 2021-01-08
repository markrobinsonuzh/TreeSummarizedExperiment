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
                       colTree = NULL, colNodeLab = NULL) {
    
    if (!is.null(rowTree)) {
        rlab <- c(rowTree$tip.label, rowTree$node.label)
        if (is.null(rowNodeLab)) {
            rowNodeLab <- rownames(x)
        }
        indR <- rowNodeLab %in% rlab
        if (sum(!indR)) {
            warning(sum(!indR), " rows are removed due to mismatch")}
        
        x <- x[indR, ]
        out <- .linkFun(tree = rowTree, 
                        sce = x,
                        nodeLab = rowNodeLab[indR],
                        onRow = TRUE)
        x@rowTree <- list(phylo = rowTree)
        x@rowLinks <- out$link
    }
    
    if (!is.null(colTree)) {
        clab <- c(colTree$tip.label, colTree$node.label)
        if (is.null(colNodeLab)) {
            colNodeLab <- colnames(x)[indC]
        }
        indC <- colNodeLab %in% clab
        if (sum(!indC)) {
            warning(sum(!indC), " cols are removed due to mismatch")}
        
        x <- x[, indC]
        out <- .linkFun(tree = colTree,
                        sce = x,
                        nodeLab = colNodeLab[indC],
                        onRow = FALSE)
        x@colTree <- list(phylo = colTree)
        x@colLinks <- out$link
    }
    
    
    return(x)
}

