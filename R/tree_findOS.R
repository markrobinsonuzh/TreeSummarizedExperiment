#' Find descendants (or offsprings)
#'
#' \code{findOS} finds descendants of a node.
#'
#' @param node An internal node. It could be the node number or the node
#'   label.
#' @param tree A phylo object.
#' @param only.leaf A logical value, TRUE or FALSE. The default is TRUE. If
#'   default, only the leaf nodes in the descendant nodes would be returned.
#' @param self.include A logical value, TRUE or FALSE. The default is FALSE. If
#'   TRUE, the node specified in \strong{node} is included and the leaf node
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
#' (tips <- findOS(tree = tinyTree, node = c(17), only.leaf = TRUE))

findOS <- function(tree,
                   node,
                   only.leaf = TRUE,
                   self.include = FALSE,
                   use.alias = FALSE) {

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }

    if (!(is.character(node) |
          is.numeric(node) |
          is.integer(node))) {
        stop("The argument (node) should be character or numeric")
    }
    # the edge matrix
    mat <- tree$edge
    matN <- matTree(tree = tree)

    if (is.character(node)) {
        numA <- transNode(tree = tree, node = node,
                          use.alias = TRUE,
                          message = FALSE)
    } else {
        numA <- node
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

    matNN <- apply(matN, 2, FUN = function(x) {
        xe <- x[!is.na(x)]
        xx <- transNode(tree = tree, node = xe,
                        use.alias = use.alias,
                        message = FALSE)
        x[!is.na(x)] <- xx
        return(x)
    })


    # descendants: tips or leaves
    tipA <- unique(setdiff(mat[, 2], mat[, 1]))

    desA <- lapply(seq_along(loc2), FUN = function(x) {
        xx <- loc2[[x]]
        y0 <- y <- matN[xx]

        # if self.include, the leaf has itself as the descendant
        if (!self.include) {
         y <- setdiff(y, numA[x])
         }

        # only leaf nodes
        if (only.leaf) {
            y <- intersect(y, tipA)
            }


        # index those kept
        ii <-  y0 %in% y
        xi <- xx[ii, , drop = FALSE]

        yy <- y0[ii]
        names(yy) <- matNN[xi]

        uy <- unique(yy)
        return(uy)
    })

    # final output (node number or label)
    names(desA) <- transNode(tree = tree, node = numA,
                             use.alias = use.alias,
                             message = FALSE)
    return(desA)
}
