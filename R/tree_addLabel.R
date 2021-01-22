#' add labels to nodes of a tree
#' 
#' \code{addLabel} label nodes of a tree (\code{phylo} object)
#' 
#' @param tree A phylo object
#' @param label A character vector to provide node labels. The label is passed to
#'   nodes that are sorted by their node number in ascending order. The default
#'   is NULL, nodes are labeled by adding a prefix \code{Node_} to their node number. 
#' @param on Chosen from "all", "leaf", "internal". If "all", all nodes are
#'   labeled; if "leaf", leaves are labeled; if "internal", internal nodes are
#'   labeled.
#' @export
#' @return a phylo object
#' @author Ruizhu Huang
#' @examples  
#' data(tinyTree)
#' library(ggtree)
#'
#' # PLOT tree
#' # The node labels are in orange texts and the node numbers are in blue
#' ggtree(tinyTree, branch.length = 'none')+
#'     geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'     geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7)
#'
#' # change labels
#' nodes <- showNode(tree = tinyTree, only.leaf = FALSE)
#' tt <- addLabel(tree = tinyTree, label = LETTERS[nodes],
#'                on = "all")
#' 
#' ggtree(tt, branch.length = 'none')+
#'     geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'     geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7)
#' 
addLabel <- function(tree, label = NULL,
                     on = c("all", "leaf", "internal")) {
    on <- match.arg(on)
    
    # edges; leaves; internal nodes
    ed <- tree$edge
    leaf <- setdiff(ed[, 2], ed[, 1])
    leaf <- sort(leaf)
    inNode <- setdiff(ed[, 1], leaf)
    inNode <- sort(inNode)
    allNode <- c(leaf, inNode)
    
    if (on == "leaf") {
        if (is.null(label)) {label <- paste0("Node_", leaf)}
        if (length(label) != length(leaf)) {
            stop("The number of labels and leaves are different.",
                 call. = FALSE)
        }
        tree <- .add_leaf_label(tree = tree, label = label)
    }
    
    if (on == "internal") {
        if (is.null(label)) {label <- paste0("Node_", inNode)}
        if (length(label) != length(inNode)) {
            stop("The number of labels and internal nodes are different.",
                 call. = FALSE)
        }
        tree <- .add_inner_label(tree = tree, label = label)
    }
   
    if (on == "all") {
        if (is.null(label)) {label <- paste0("Node_", c(leaf, inNode))}
        if (length(label) != length(allNode)) {
            stop("The number of labels and nodes are different.",
                 call. = FALSE)
        }
        leafLab <- label[seq_along(leaf)]
        inLab <- setdiff(label, leafLab)
        tree <- .add_leaf_label(tree = tree, label = label[seq_along(leaf)])
        tree <- .add_inner_label(tree = tree, label = inLab)
    }
    
    lab <- c(tree$tip.label, tree$node.label)
    if (anyDuplicated(lab)) {
        warning("The tree has duplicated labels")
    }
    return(tree)
}

.add_leaf_label <- function(tree, label) {
    tree$tip.label <- label
    return(tree)
}
.add_inner_label <- function(tree, label) {
    tree$node.label <- label
    return(tree)
}

