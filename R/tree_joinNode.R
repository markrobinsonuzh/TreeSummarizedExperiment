#' find the optimal nodes to short result.
#'
#' \code{joinNode} is to use as few as possible nodes to represent the provided
#' nodes so that descendant leaves covered by the input nodes and output nodes
#' are exactly the same.
#'
#' @param tree A tree (phylo object)
#' @param node A vector of node numbers or node labels
#' @param use.alias A logical value, TRUE or FALSE. The default is FALSE, and
#'   the node label would be used to name the output; otherwise, the alias of
#'   node label would be used to name the output. The alias of node label is
#'   created by adding a prefix \code{"alias_"} to the node number.
#' @export
#' @return A vector of nodes. The numeric value is the node number, and the
#'   vector name is the corresponding node label. If a node has no label, it
#'   would have NA as name when \code{use.alias = FALSE}, and have the alias of
#'   node label as name when \code{use.alias = TRUE}.
#' @author Ruizhu Huang
#' @examples
#'
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
#' ## find the node shared by provided node labels
#' joinNode(node = c('t4','t9'), tree = tinyTree)
#' joinNode(node = c('t4','t9'), tree = tinyTree)
#' joinNode(node = c('t10','Node_18', 't8'), 
#'          tree = tinyTree,
#'          use.alias = FALSE)
#' joinNode(node = c('t10','Node_18', 't8'), 
#'          tree = tinyTree,
#'          use.alias = TRUE)
#'
#' ## find the node shared by provided node numbers
#' joinNode(node = c(2, 3), tree = tinyTree)
#' joinNode(node = c(2, 3, 16), tree = tinyTree)
#'
joinNode <- function(tree, node, use.alias = FALSE) {

    if (!inherits(tree, "phylo")) {
        stop("tree is not a phylo object.")
    }

    if (!is.atomic(node)) {
        stop("node is a vector")
    }

    # transfer node label to node number
    if (is.character(node)) {
        node <- convertNode(tree, node = node,
                          message = FALSE)
    } else {
        node <- node
    }

    # path matrix
    mat <- matTree(tree)

    if (!all(node %in% as.vector(mat))) {
        stop("Some nodes could not be found in the tree")
    }

    # select paths which include the input nodes
    ind <- apply(mat, 1, FUN = function(x) {
        any(x %in% node)
    })
    # select nodes which only exist in the selected paths
    selN <- setdiff(as.vector(mat[ind, ]), as.vector(mat[!ind, ]))
    selN <- selN[!is.na(selN)]
    # remove nodes which are descendants of any others
    matI <- mat[ind,, drop = FALSE]
    selF <- apply(matI, MARGIN = 1,
                FUN = function(x){
                    y <- x %in% selN
                    x[sum(y)]
                })
    sNode <- unique(selF)
    out <- sNode[!is.na(sNode)]

    # final output
    names(out) <- convertNode(tree = tree, node = out,
                            use.alias = use.alias,
                            message = FALSE)
    return(out)
}
