#' To print out the node labels
#'
#' \code{nodeLabel} is to print out the node labels of a \code{phylo} tree.
#'
#' @param tree A phylo object.
#' @param type A character value choose from \code{leaf}, \code{all}, and
#'   \code{internal}. If \code{leaf}, the output is a data frame including only
#'   leaf nodes; if \code{internal}, the output is a data frame including only
#'   internal nodes; if \code{all}, the output is a data frame including all
#'   nodes.
#'
#' @export
#' @author Ruizhu HUANG
#' @return a data frame
#' @examples
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
#' (pn1 <- printNode(tinyTree, type = "leaf"))
#' (pn2 <- printNode(tinyTree, type = "internal"))
#' (pn3 <- printNode(tinyTree, type = "all"))
#'
printNode <- function(tree, type = "leaf"){

    # edges
    edge <- tree$edge

    # the alias label
    labAlias <- tree$alias.label

    if (is.null(labAlias)) {
        message("The alias label can't be found;
                 Prefix 'alias_' to the node number as the alias label.")

    }
    # leaves
    leaf <- setdiff(edge[, 2], edge[, 1])
    leaf <- sort(leaf)
    labL <- tree$tip.label

    outL <- data.frame(nodeLab = labL,
                       nodeLab_alias = paste0("alias", leaf),
                       nodeNum = leaf,
                       isLeaf = TRUE,
                       stringsAsFactors = FALSE)


    # all nodes
    nodeA <- unique(as.vector(edge))
    nodeA <- sort(nodeA)


    # internal nodes
    nodeI <- setdiff(nodeA, leaf)
    labI <- tree$node.label
    outI <- data.frame(nodeLab_alias = paste0("alias_", nodeI),
                       nodeNum = nodeI,
                       isLeaf = FALSE)
    if (length(labI) == length(nodeI)) {
            outI$nodeLab <- labI
        } else{
            if (length(labI)) {
                warning("The node.label is less than the internal nodes")
                outI$nodeLab <- rep_len(labI, length(nodeI))
            } else {
                outI$nodeLab <- NA
            }
        }
    outI <- outI[, c("nodeLab", "nodeLab_alias", "nodeNum", "isLeaf")]

    if (type == "leaf") { out <- outL }
    if (type == "internal") { out <- outI }
    if (type == "all") { out <- rbind(outL, outI)}

    return(out)
}
