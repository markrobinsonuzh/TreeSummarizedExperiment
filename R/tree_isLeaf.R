#' To test whether the specified nodes are leaf nodes
#'
#' \code{isLeaf} is to test wheter some specified nodes are leaf nodes of a
#' tree.
#'
#' @param tree A phylo object.
#' @param node A numeric or character vector. Node labels or node numbers.
#'
#' @export
#' @author Ruizhu HUANG
#' @return a logical vector with the same length as the input \code{node}.
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
#' isLeaf(tree = tinyTree, node = c(5, 4, 18))
#' isLeaf(tree = tinyTree, node = c("t4", "t9", "Node_18" ))
#'
isLeaf <- function(tree, node){

    edge <- tree$edge
    leaf <- setdiff(edge[, 2], edge[, 1])

    # the input node
    if (is.character(node)) {
        node <- convertNode(tree, node = node,
                           message = FALSE)
    } else {
        if (is.numeric(node)) {
            node <- node
        } else {
            stop("The input node should be a character or numeric vector. \n")
        }
    }

    # check whether the specified node exists on the tree
    if (any(!node %in% as.vector(edge))) {
        stop("Node", node[!node %in% as.vector(edge)],
             " can't be matched to any node of the tree. \n")
    }


    # is it a leaf node
    isLeaf <- node %in% leaf
    return(isLeaf)
}
