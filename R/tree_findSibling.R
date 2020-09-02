#' find the sibling node
#'
#' \code{findSibling} is to find the sibling node of an \code{node} node.
#'
#' @param tree A phylo object.
#' @param node A numeric or character vector. Node labels or node numbers.
#' @param use.alias A logical value, TRUE or FALSE. The default is FALSE, and
#'   the original node label would be used to name the output; otherwise, the
#'   alias of node label would be used to name the output. The alias of node
#'   label is created by adding a prefix \code{"alias_"} to the node number.
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

    if (is.character(node)) {
        node <- convertNode(tree = tree, node = node, message = FALSE)
    }
    matT <- matTree(tree = tree)

    # find the parent node of the input node
    pN <- findAncestor(tree = tree, node = node, level = 1)

    # find the siblings
    loc1 <- lapply(pN, FUN = function(x) {
        which(matT == x, arr.ind = TRUE)
    })

    loc2 <- lapply(loc1, FUN = function(x) {
        x1 <- x[, "col"]
        x2 <- x1 - 1

        mc <- cbind(row = x[, "row"],
                    col = x2)
    })

    sib <- lapply(seq_along(loc2), FUN = function(x){
        xx <- matT[loc2[[x]]]
        setdiff(xx, node[x])

        })
    out <- unlist(sib)

    names(out) <- convertNode(tree = tree, node = out,
                            use.alias = use.alias,
                            message = FALSE)
    return(out)
}
