#' Find the ancestors of specified nodes
#'
#' \code{findAncestor} finds the ancestor in the nth generation above
#' specified nodes.
#'
#' @param tree A phylo object
#' @param node A vector of node numbers or node labels
#' @param level A vector of numbers to define nth generation before the
#' specified nodes
#' @param use.alias A logical value, TRUE or FALSE. The default is FALSE, and
#'   the node label would be used to name the output; otherwise, the alias of
#'   node label would be used to name the output. The alias of node label is
#'   created by adding a prefix \code{"Node_"} to the node number if the node is
#'   an internal node or adding a prefix \code{"Leaf_"} if the node is a leaf
#'   node.
#' @export
#' @return A vector of nodes. The numeric value is the node number, and the
#'   vector name is the corresponding node label. If a node has no label, it
#'   would have NA as name when \code{use.alias = FALSE}, and have the alias of
#'   node label as name when \code{use.alias = TRUE}.
#' @author Ruizhu Huang
#'
#' @examples
#' library(ggtree)
#' data(tinyTree)
#' ggtree(tinyTree, branch.length = 'none') +
#'  geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'  geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7)
#'
#'  findAncestor(tree = tinyTree, node = c(18, 13), level = 1)

findAncestor <- function(tree, node, level,
                         use.alias = FALSE) {

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }

    # convert a tree to a matrix
    # each row is a path connecting the root and a leaf
    treeMat <- matTree(tree)


    if (is.character(node)) {
        node <- transNode(tree = tree, node = node,
                            use.alias = TRUE,
                            message = FALSE)
    }

    if (length(level) == 1) {
        level <- rep(level, length(node))
    } else {
        if (length(level) == length(node)) {
            level <- level
        } else {
            stop("the length of level is not equal to the length of node")
        }
    }

    selNod <- lapply(seq_along(node), FUN = function(x) {
        # the node
        nod.x <- node[x]

        # where is the node
        ind <- which(treeMat == nod.x, arr.ind = TRUE)

        # the level
        level.x <- level[x]
        ind.x <- ind
        ind.x[, "col"] <- ind[, "col"] + level.x

        if (any (ind.x[, "col"] > ncol(treeMat))) {
            stop("Exceed the root; try a lower level.")
        }

        vv <- treeMat[ind.x]
        uv <- unique(as.vector(vv))

        if (length(uv) > 1) {
            stop("More than one node are found.")
        }
        return(uv)
    })

    out <- unlist(selNod)

    # return a vector of the found node (the node number of the node)
    # name the vector with the node label
    names(out) <- transNode(tree = tree, node = out,
                            use.alias = use.alias,
                            message = FALSE)
    return(out)
}
