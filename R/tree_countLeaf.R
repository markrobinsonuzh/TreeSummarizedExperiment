#' count the number of leaf nodes
#'
#' \code{countLeaf} calculates the number of leaves on a \code{phylo} tree.
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


countLeaf <- function(tree) {

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object. \n")
    }


    # node number & tip number
    mat <- tree$edge
    tip <- setdiff(mat[, 2], mat[, 1])
    utip <- unique(tip)

    len <- length(tip)

    return(len)

    }

