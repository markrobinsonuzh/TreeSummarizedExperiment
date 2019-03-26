#' find the sibling node
#'
#' \code{findSibling} is to find the sibling node of an \code{node} node.
#'
#' @param tree A phylo object.
#' @param node A numeric or character vector. Node labels or node numbers.
#' @param use.alias A logical value, TRUE or FALSE. The default is FALSE, and
#'   the original node label would be used to name the output; otherwise, the
#'   alias of node label would be used to name the output. The alias of node
#'   label is created by adding a prefix \code{"Node_"} to the node number if
#'   the node is an internal node or adding a prefix \code{"Leaf_"} if the node
#'   is a leaf node.
#' @export
#' @return A vector of nodes. The numeric value is the node number, and the
#'   vector name is the corresponding node label. If a node has no label, it
#'   would have NA as name when \code{use.alias = FALSE}, and have the alias of
#'   node label as name when \code{use.alias = TRUE}.
#'
#' @examples
#' library(ggtree)
#' data(tinyTree)
#' ggtree(tinyTree, branch.length = 'none') +
#'     geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'     geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7)
#'
#'
#'  findSibling(tree = tinyTree, node = 17)
#'  findSibling(tree = tinyTree, node = c(13, 17))
findSibling <- function(tree, node, use.alias = FALSE){

    # find descendant leaves of the input node
    inT <- findOS(tree = tree, node = node,
               only.leaf = TRUE)
    # find the parent node of the input node
    pN <- findAncestor(tree = tree, node = node, level = 1)

    # Leaves not included in input node
    exT <- lapply(seq_along(pN), FUN = function(x) {
        aT <- findOS(tree = tree, node = pN[x], only.leaf = TRUE)[[1]]
        setdiff(aT, inT[[x]])
    })

    # replace leaves with their ancestor branch node
    fT <- lapply(exT, FUN = function(x) {
        signalNode(tree = tree, node = x,
                   use.alias = use.alias)} )
    out <- unlist(fT)

    names(out) <- transNode(tree = tree, node = out,
                            use.alias = use.alias,
                            message = FALSE)
    return(out)
}
