#' Transfer between node number and node label
#'
#' \code{transNode} does the transformation between the number and the label of
#' a node on a tree
#'
#' @param tree A phylo object
#' @param node A character or numeric vector representing tree node label(s) or
#'   tree node number(s)
#' @param use.alias A logical value, TRUE or FALSE. If the input \code{node} is
#'   character and \code{use.alias = TRUE}, it tells that the input label is the
#'   alias of node labels; otherwise, it tells that the input label is the
#'   original node labels; If the input \code{node} is numeric and
#'   \code{use.alias = TRUE}, the alias of node labels are returned; otherwise,
#'   the original node labels are returned.
#' @param message A logical value, TRUE or FALSE. The default is FALSE. If TRUE,
#'   message will show when a tree have duplicated labels for some internal
#'   nodes.
#'
#' @export
#' @return a vector
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
#' #check whether the node number and node label are matched
#' transNode(tinyTree, node = c(11, 2, 4, 15), use.alias = FALSE)
#' transNode(tinyTree, node = c(11, 2, 4, 15), use.alias = TRUE)
#'
#' transNode(tree = tinyTree, node = c("Node_16", "Node_11"))
#' transNode(tree = tinyTree, node = c("alias_16", "alias_11"),
#' use.alias = TRUE)

transNode <- function(tree, node, use.alias = FALSE,
                      message = FALSE) {

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object. \n")
    }


    if (is.factor(node)) {
        stop("factor detected; The node label is required to be character or
             numeric.")
    }
    # node number & tip number
    mat <- tree$edge
    nodI <- sort(unique(mat[, 1]))
    tip <- sort(setdiff(mat[, 2], mat[, 1]))
    nodeA <- sort(c(tip, nodI))
    if (is.null(tree$alias.label)) {
        tree$alias.label <- paste0("alias_", nodeA)
    }

    # if node labels are given, check whether the length could match with the
    # length of internal nodes.
    if (!is.null(tree$node.label)) {
        if (length(tree$node.label) != length(nodI)) {
            stop("The length of internal node label isn't equal to
                 the length of the internal nodes. \n")
        }
    }

    # node labels
    nodeLab <- c(tree$tip.label, tree$node.label)
    nodeLab_alias <- tree$alias.label
    if (message) {
        if (any(duplicated(nodeLab))) {
            cat("There are more than one nodes using a same label or
                without any label.\n")
        }
    }

    # check whether the input node number exists in the provided tree
    if (is.numeric(node)) {
        if (!all(node %in% nodeA)) {
            stop("The node number ", node[!node %in% nodeA],
                 " can't be found in the ",
                 deparse(substitute(tree)), "\n")
        }
    }
    # check whether the input label exists in the provided tree
    # (allow nodeLab_alias)
    inLab <- node %in% nodeLab
    inAlias <- node %in% nodeLab_alias

    if (is.character(node)) {
        if (use.alias) {
            if (!all(inAlias)) {
                mlen <- ifelse(sum(!inAlias) > 4, 4, sum(!inAlias))
                mis <- node[!inAlias][seq_len(mlen)]
                stop("Can't find ", paste(mis, collapse = " "))
            }
        } else {
            if (!all(inLab)) {
                mlen <- ifelse(sum(!inLab) > 4, 4, sum(!inLab))
                mis <- node[!inLab][seq_len(mlen)]
                if (all(startsWith(mis, "alias"))) {
                    message("You might want to try 'use.alias = TRUE'")
                }
                stop("Can't find ", paste(mis, collapse = " "))
            }
        }
     }

    # =============== Translation ======================
    # translate from the label to the number
    if (is.character(node)) {
        if (use.alias) {
            names(nodeA) <- nodeLab_alias
            final <- nodeA[node]
        } else {
            names(nodeA) <- nodeLab
            final <- nodeA[node]
        }
    }

    # transfer from the number to the label
    if (is.numeric(node)) {
        if (use.alias) {
            sel <- match(node, nodeA)
            final <- nodeLab_alias[sel]
        } else {
            sel <- match(node, nodeA)
            final <- nodeLab[sel]
        }}

    # output
    return(final)

}

