#' Find the share node
#'
#' \code{shareNode} is to find the node where the specified nodes first meet.
#'
#' @param tree A  phylo object.
#' @param node A vector of node numbers or node labels.
#' @param use.alias A logical value, TRUE or FALSE. The default is FALSE, and
#'   the node label would be used to name the output; otherwise, the alias of
#'   node label would be used to name the output. The alias of node label is
#'   created by adding a prefix \code{"Node_"} to the node number if the node is
#'   an internal node or adding a prefix \code{"Leaf_"} if the node is a leaf
#'   node.
#'
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
#' shareNode(node = c('t4','t9'), tree = tinyTree,
#'           use.alias = FALSE)
#'
#' shareNode(node = c('t10','Node_17'), tree = tinyTree,
#'           use.alias = FALSE)
#'
#' ## find the node shared by provided node numbers
#' shareNode(node = c(2, 3), tree = tinyTree)

shareNode <- function(tree, node,
                      use.alias = FALSE) {

    if (!inherits(tree, "phylo")) {
        stop("tree is not a phylo object.")
    }

    if (!is.atomic(node)) {
        stop("node is a vector")
    }

    # transfer node label to node number
    if (is.character(node)) {
        node <- transNode(tree, node = node,
                          message = FALSE)
    } else {
        node <- node
    }

    # path matrix
    mat <- matTree(tree)
    ind <- apply(mat, 1, FUN = function(x) {
        any(x %in% node)
    })
    matN <- mat[ind, ]
    path <- lapply(seq_len(nrow(matN)), FUN = function(x) {
        xx <- matN[x, ]
        xx[!is.na(xx)]
    })
    # ancestors
    loc <- Reduce(intersect, path)

    # the ancestor on the lowest level (the root has the highest level)
    vec <- as.vector(mat)
    count <- table(vec, useNA = "no")[loc]
    if (length(count) == 1) {
        af <- as.data.frame(count, stringsAsFactors = FALSE)
        df <- data.frame(vec = as.numeric(rownames(af)),
                         Freq = af$count)
    } else {
        df <- as.data.frame(count, stringsAsFactors = FALSE)
    }


    # select the node with the lowest frequency.  closest to the leaf level.
    out <- as.numeric(df[df$Freq == min(df$Freq), "vec"])

    # final output
    names(out) <- transNode(tree = tree, node = out,
                            use.alias = use.alias,
                            message = FALSE)
    return(out)
   }
