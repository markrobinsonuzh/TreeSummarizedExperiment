#' Find nodes on the tree
#'
#' \code{showNode} is to get nodes from the tree.
#'
#' @param tree A  phylo object.
#' @param only.leaf A logical value, TRUE or FALSE. The default is FALSE, all
#'   nodes are output; otherwise, leaves are output
#' @param use.alias A logical value, TRUE or FALSE. The default is FALSE, and
#'   the node label would be used to name the output; otherwise, the alias of
#'   node label would be used to name the output. The alias of node label is
#'   created by adding a prefix \code{"alias_"} to the node number.
#' @export
#' @importFrom methods is
#' @return A vector of nodes. The numeric value is the node number, and the
#'   vector name is the corresponding node label. If a node has no label, it
#'   would have NA as name when \code{use.alias = FALSE}, and have the alias of
#'   node label as name when \code{use.alias = TRUE}.
#' @author Ruizhu Huang
#'
#' @examples
#' library(ggtree)
#' data(tinyTree)
#'
#' # PLOT tree
#' ggtree(tinyTree, branch.length = 'none') +
#'     geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'     geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7)
#'
#'
#' ## find the node shared by provided node labels
#' showNode(tree = tinyTree, only.leaf = TRUE,
#'           use.alias = FALSE)
#'
#' showNode(tree = tinyTree, only.leaf = FALSE,
#'           use.alias = FALSE)

showNode <- function(tree, only.leaf = FALSE,
                      use.alias = FALSE) {
    
    if (!is(tree, "phylo")) {
        stop("tree is not a phylo object.")
    }
    
    # edges
    ed <- tree$edge
    
    # leaves
    leaf <- setdiff(ed[, 2], ed[, 1])
    leaf <- unique(leaf)
    leaf <- sort(leaf)
    
    # nodes
    all <- as.vector(ed)
    all <- unique(all)
    all <- sort(all)
    
    # output
    if (only.leaf) {
        out <- leaf
    } else {
        out <- all
    }
    
    # final output
    names(out) <- transNode(tree = tree, node = out,
                            use.alias = use.alias,
                            message = FALSE)
    return(out)
}
