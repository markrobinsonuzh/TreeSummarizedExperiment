#' track the nodes of a phylo tree
#'
#' \code{trackNode} track nodes of a phylo tree by adding the alias labels to
#' them
#'
#' @param tree A phylo object
#'
#' @export
#' @return a phylo object
#' @author Ruizhu Huang
#'
#' @examples
#' library(ggtree)
#'
#' data(tinyTree)
#'
#' ggtree(tinyTree, branch.length = 'none') +
#'     geom_text2(aes(label = label), hjust = -0.3) +
#'     geom_text2(aes(label = node), vjust = -0.8,
#'                hjust = -0.3, color = 'blue')
#'
#' #check whether the node number and node label are matched
#' trackTree <- trackNode(tinyTree)
#' ggtree(trackTree, branch.length = 'none') +
#'     geom_text2(aes(label = label), hjust = -0.3) +
#'     geom_text2(aes(label = node), vjust = -0.8,
#'                hjust = -0.3, color = 'blue')

trackNode <- function(tree) {
    
    # ===========check inputs =================
    if (!is(tree, "phylo")) {
        stop("tree should be a phylo object")
    }
    
    # The node number of all nodes
    ed <- tree$edge
    leaf <- setdiff(ed[, 2], ed[, 1])
    leaf <- sort(leaf)
    inNode <- unique(sort(ed[, 1]))


    # prefix node number with 'alias_' to create the alias label
    tree$tip.label <- paste0("alias_", leaf)
    tree$node.label <- paste0("alias_", inNode)
    return(tree)
}
