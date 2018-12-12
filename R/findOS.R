#' Find descendants (or offsprings)
#'
#' \code{findOS} finds descendants of a node.
#'
#' @param ancestor An internal node. It could be the node number or the node
#'   label.
#' @param tree A phylo object.
#' @param only.Tip A logical value, TRUE or FALSE. The default is TRUE. If
#'   default, only the leaf nodes in the descendant nodes would be returned.
#' @param self.include A logical value, TRUE or FALSE. The default is TRUE. If
#'   default, the node specified in \strong{ancestor} is included. The leaf node
#'   itself is returned as its descendant.
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
#' data(tinyTree)
#'
#' library(ggtree)
#' ggtree(tinyTree) +
#' geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7) +
#' geom_hilight(node = 17, fill = 'steelblue', alpha = 0.5) +
#' geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7)
#'
#' (tips <- findOS(tree = tinyTree, ancestor = c(17), only.Tip = TRUE))

findOS <- function(tree,
                   ancestor,
                   only.Tip = TRUE,
                   self.include = TRUE,
                   use.alias = FALSE) {

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }

    if (!(is.character(ancestor) |
          is.numeric(ancestor) |
          is.integer(ancestor))) {
        stop("ancestor should be character or numeric")
    }
    # the edge matrix
    mat <- tree$edge
    matN <- matTree(tree = tree)

    if (is.character(ancestor)) {
        numA <- transNode(tree = tree, input = ancestor,
                          use.alias = TRUE,
                          message = FALSE)
    } else {
        numA <- ancestor
        isOut <- !numA %in% mat
        if (any(isOut)) {
            stop("Node ", numA,
                 " can't be found in the ",
                 deparse(substitute(tree)), "\n")
            }

    }

    # find all descendants
    loc1 <- lapply(numA, FUN = function(x) {
        xi <- which(matN == x, arr.ind = TRUE)
        return(xi)
    })

    loc2 <- lapply(loc1, FUN = function(x) {
        x1 <- lapply(x[, "col"], seq)
        x1 <- unlist(x1)
        x2 <- rep(x[, "row"], x[, "col"])
        cbind(row = x2, col = x1)

    })

    desA <- lapply(loc2, FUN = function(x) {
        matN[x]
    })


    # descendants: tips or leaves
    tipA <- unique(setdiff(mat[, 2], mat[, 1]))

    # if self.include, the leaf has itself as the descendant
    if (!self.include) {
        desA <- lapply(seq_along(numA), FUN = function(x) {
            setdiff(desA[[x]], numA[x])
        })
    }

    # if only.Tip, the descendant leaves are return and the descendant internal
    # nodes are not
    if (only.Tip) {
        desA <- lapply(desA, FUN = function(x) {
            intersect(x, tipA)
        })
    }


    # name the node number with the node label
    desA <- lapply(desA, FUN = function(x) {
        xx <- transNode(tree = tree, input = x,
                        use.alias = use.alias,
                        message = FALSE)
        names(x) <- xx
        return(x)
    })


    # final output (node number or label)
    names(desA) <- transNode(tree = tree, input = numA,
                            use.alias = use.alias,
                            message = FALSE)
    return(desA)
}
