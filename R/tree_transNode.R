#' Transfer between node number and node label
#'
#' \code{transNode} does the transformation between the number and the label of
#' a node on a tree
#'
#' @param tree A phylo object
#' @param node A character or numeric vector representing tree node label(s) or
#'   tree node number(s)
#' @param use.alias A logical value, TRUE or FALSE. This is an optional argument
#'   that only requried when the input \code{node} is a numeric vector. The
#'   default is FALSE, and the node label would be returned; otherwise, the
#'   alias of node label would be output. The alias of node label is created by
#'   adding a prefix \code{"alias_"} to the node number.
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
#' transNode(tinyTree, node = c(11, 2, 4, 15))
#'
#' transNode(tree = tinyTree, node = c("Node_16", "Node_11"))
#' transNode(tree = tinyTree, node = c("alias_16", "alias_11"))

transNode <- function(tree, node, use.alias = FALSE,
                      message = FALSE) {

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object. \n")
    }

    if (is.factor(node)) {
        node <- as.character(node)
        warning("node: factor is converted to character")
        warning("The input 'node' is considered as labels of nodes.")
    }
    # node number & tip number
    mat <- tree$edge
    nodI <- sort(unique(mat[, 1]))
    tip <- sort(setdiff(mat[, 2], mat[, 1]))
    nodeA <- c(tip, nodI)

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
    nodeLab_alias <- paste("alias_", c(tip, nodI), sep = "")
    
    if (any(duplicated(nodeLab))) {
        warnings("Multiple nodes use the same label or
                have no label.\n")
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
    inLab <- all(node %in% nodeLab)
    inAlias <- all(node %in% nodeLab_alias)
    if (is.character(node)) {
        if (!any(inLab, inAlias)) {
            wrongNode <- setdiff(node, nodeLab)
            if (any(startsWith(wrongNode, "alias_"))) {
                message("Not allowed to mix using node labels and alias labels")
            }
            stop(length(wrongNode), " nodes mismatch with the tree: ",
                    head(wrongNode), " ...")
            

        }
    }

    # =============== Transformation ======================
    # transfer from the label to the number
    if (is.character(node)) {
        if (inLab) {
            names(nodeA) <- nodeLab
            final <- nodeA[node]
        } else {
            names(nodeA) <- nodeLab_alias
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
