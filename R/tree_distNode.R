#' Calculate the distance between any two nodes on the tree
#'
#' \code{distNode} is to calculate the distance between any two nodes on
#' a \code{phylo} tree
#'
#' @param tree A phylo object.
#' @param node A numeric or character vector of length two.
#' @export
#' @return A numeric value.
#'
#' @examples
#' library(ggtree)
#' data(tinyTree)
#' ggtree(tinyTree) +
#'     geom_text2(aes(label = node), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'     geom_text2(aes(label = branch.length), color = "darkblue",
#'                vjust = 0.7)
#'
#'
#' distNode(tree = tinyTree, node = c(10, 11))
#' distNode(tree = tinyTree, node = c(12, 13))
#' distNode(tree = tinyTree, node = c(13, 15))
#' distNode(tree = tinyTree, node = c(12, 14))

distNode <- function(tree, node) {
    
    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object. \n")
    }
    
    if (is.character(node)) {
        node <- convertNode(tree = tree, node = node)
    }

    mat <- matTree(tree)
    ed <- tree$edge

    oo <- lapply(node, FUN = function(x) {
        xx <- which(mat == x, arr.ind = TRUE)
        return(xx)
    })

    samePath <- intersect(oo[[1]][, "row"], oo[[2]][, "row"])
    if (length(samePath)) {
        dd <- .samePath(tree = tree, node = node)
    } else {
        dd <- .diffPath(tree = tree, node = node)
    }

    return(dd)

    }

.samePath <- function(tree, node) {
    # edges
    ed <- tree$edge

    # edge.length
    len <- tree$edge.length

    # the input nodes are not directly connected but in the same path
    mt <- matTree(tree = tree)
    loc <- lapply(node, FUN = function(x) {
        xx <- which(mt == x, arr.ind = TRUE)
    })

    shareRow <- intersect(loc[[1]][, "row"], loc[[2]][, "row"])
    ush <- unique(shareRow)[1]

    loc <- do.call(rbind, loc)
    tcol <- loc[loc[, "row"] == ush, "col"]
    sel <- seq(from = max(tcol), to = min(tcol), by = -1)
    pn <- mt[ush, sel[-length(sel)]]
    cn <- mt[ush, sel[-1]]
    pir <- cbind(pn, cn)
    fi <- match(data.frame(t(pir)), data.frame(t(ed)))
    dd <- sum(len[fi])

    return(dd)
}

.diffPath <- function(tree, node) {
    # the share node
    sNode <- shareNode(tree = tree, node = node)


    d1 <- .samePath(tree = tree, node = c(node[1], sNode))
    d2 <- .samePath(tree = tree, node = c(node[2], sNode))
    dd <- d1 + d2

    return(dd)
}


