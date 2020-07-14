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
    num <- transNode(tree = newTree, node = newLab)
    lab <- transNode(tree = newTree, node = num)
    alias <- transNode(tree = newTree,
                       node = num, 
                       use.alias = TRUE)
    ind <- isLeaf(tree = newTree, node = num)
    LinkDataFrame(nodeLab = lab, nodeLab_alias = alias, 
                  nodeNum = num, isLeaf = ind)
}
# updateTree <- function(x, tree, on_row = TRUE,
#                        match_node = NULL, message = TRUE) {
#     
#     # the old tree
#     if (on_row) {
#         oldTree <- rowTree(x)
#         oldLink <- DataFrame(rowLinks(x))
#     } else {
#         oldTree <- colTree(x)
#         oldLink <- DataFrame(colLinks(x))
#     }
#    
#     # nodes: in the old tree but not in the new tree
#     oldLabel <- c(oldTree$tip.label, oldTree$node.label)
#     newLabel <- c(tree$tip.label, tree$node.label)
#     diffLabel <- setdiff(oldLabel, newLabel)
#     
#     if (!is.null(match_node)) {
#         # each column stores the node number
#         match_class <- unique(apply(match_node, 2, class))
#         if (length(match_class) > 1) {
#             stop("match_node should have the same class in two columns")
#         }
#         if (!is.numeric(match_class)) {
#             match_node <- mapply(transNode, list(new = tree, old = oldTree), 
#                                  match_node)
#             match_node <- match_node[!duplicated(match_node), , drop = FALSE]
#         }
#         if (anyDuplicated(match_node[, "new"])) {
#             Stop("Multiple new nodes are mapped to an old node")
#         }
#         matchLabel <- transNode(tree = oldTree, node = match_node[, "old"])
#         diffLabel <- setdiff(diffLabel, matchLabel)
#         
#         # remove rows/cols that don't belong to the provided nodes
#         ind <- !(oldLink$nodeNum %in% match_node[, "old"])
#         
#     } else {
#          # remove rows/cols mapped to these missing nodes
#         ind <- oldLink$nodeLab %in% diffLabel
#     }
#     if (message) {
#         message(length(diffLabel), 
#         " old nodes (with labels) can't be found in the new tree")
#     }
#     
#     
#     
#     if (on_row) {
#         xx <- x[!ind, ]   
#         rn <- rownames(xx)
#         newLink <- rowLinks(xx)
#         if (message) {
#             message(sum(ind), 
#                     " rows mismatch with nodes and are removed")}
#     } else {
#         xx <- x[, !ind]
#         rn <- colnames(xx)
#         newLink <- colLinks(xx)
#         if (message) {
#             message(sum(ind), 
#                     " cols mismatch with nodes and are removed")}
#     }
#     newLink <- DataFrame(newLink)
#     
#     # update links
#     if (!is.null(match_node)) {
#         match_alias <- mapply(transNode, list(new = tree, old = oldTree), 
#                               data.frame(match_node[, c("new", "old")]),
#                               use.alias = TRUE)
#         pair_alias <- setNames(match_alias[, "new"], match_alias[, "old"])
#         newLink$nodeLab_alias <- pair_alias[newLink$nodeLab_alias] 
#         newLink$nodeNum <- transNode(tree = tree, 
#                                      node = newLink$nodeLab_alias,
#                                      message = FALSE)
#     } else {
#         newLink$nodeNum <- transNode(tree = tree, node = newLink$nodeLab,
#                                      message = FALSE)
#         newLink$nodeLab_alias <- transNode(tree = tree, 
#                                            node = newLink$nodeNum,
#                                            use.alias = TRUE, 
#                                            message = FALSE)
#     }
#     newLink$isLeaf <- isLeaf(tree = tree, node = newLink$nodeNum)
#     newLink <- as(newLink, "LinkDataFrame")
#     
#     ## update the tree and links
#     ## update row/col names if they are alias name
#     is_alias <- startsWith(rn, "alias_")
#     
#     if (on_row) {
#         newTse <- BiocGenerics:::replaceSlots(xx,
#                                               rowLinks = newLink,
#                                               rowTree = list(phylo = tree))
#         if (any(is_alias)) {
#             rownames(newTse) <- newLink$nodeLab_alias
#         }
#     } else {
#         newTse <- BiocGenerics:::replaceSlots(xx,
#                                               colLinks = newLink,
#                                               colTree = list(phylo = tree))
#         if (any(is_alias)) {
#             colnames(newTse) <- newLink$nodeLab_alias
#         }
#     }
#     
#     return(newTse)
# }
