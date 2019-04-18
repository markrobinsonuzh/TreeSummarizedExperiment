#' count the number of nodes
#'
#' \code{countNode} calculates the number of nodes on a \code{phylo} tree.
#'
#' @param tree A phylo object
#'
#' @export
#' @return a numeric value
#' @author Ruizhu Huang
#'
#' @examples
#' library(ggtree)
#'
#' data(tinyTree)
#'
#' ggtree(tinyTree, branch.length = 'none') +
#' geom_text2(aes(label = label), hjust = -0.3) +
#' geom_text2(aes(label = node), vjust = -0.8,
#' hjust = -0.3, color = 'blue')
#'
#'
#' (n <- countLeaf(tinyTree))
#'


countNode <- function(tree) {

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object. \n")
    }


    # all nodes
    mat <- tree$edge
    node <- as.vector(mat)
    unode <- unique(node)

    len <- length(unode)

    return(len)

}

