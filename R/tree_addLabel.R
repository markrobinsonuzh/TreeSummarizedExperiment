#' add labels to nodes of a tree
#' 
#' \code{addLabel} adds labels to the node of a tree (\code{phylo} object)
#' 
#' @param tree A phylo object
#' @param label A character vector as the label of tree. The label is passed to
#'   nodes that are sorted by their node number in ascending order. The default
#'   is NULL, and nodes are labeled by using their node numbers (convert node
#'   numbers from numeric values to characters)
#'   
#' @export
#' @return a phylo object
#' @author Ruizhu Huang
#' @examples  
#' data(tinyTree)
#' library(ggtree)
#'
#' # PLOT tree
#' # The node labels are in orange texts and the node numbers are in blue
#' ggtree(tinyTree,branch.length = 'none')+
#'     geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'     geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7)
#'
#' # change labels
#' nodes <- showNode(tree = tinyTree, only.leaf = FALSE)
#' tt <- addLabel(tree = tinyTree, label = LETTERS[nodes])
#' 
#' ggtree(tt, branch.length = 'none')+
#'     geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'     geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7)
#' 
addLabel <- function(tree, label = NULL) {
    ed <- tree$edge
    leaf <- setdiff(ed[, 2], ed[, 1])
    leaf <- sort(leaf)
    inNode <- setdiff(ed[, 1], leaf)
    inNode <- sort(inNode)
    
    if (is.null(label)) {
        tree$tip.label <- as.character(leaf)
        tree$node.label <- as.character(inNode)
    } else {
        tree$tip.label <- label[seq_along(leaf)]
        tree$node.label <- setdiff(label, tree$tip.label)
    }
    
    return(tree)
}
