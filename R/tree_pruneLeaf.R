#' remove branches of a phylo tree
#'
#' \code{pruneTree} is to remove branches that the specified leaves are in from
#' a \code{phylo} tree
#'
#' @param tree A phylo object.
#' @param rmLeaf A numeric or character vector. The labels or numbers of leaves
#'   to be removed.
#' @param mergeSingle A logical value, TRUE or FALSE. If TRUE, the internal node
#'   that has only one child will be omitted. For example, a tree has structure
#'   as A-[B,C-D] (A has two children B, and C. C has a child D). If
#'   \code{mergeSingle = TRUE}, the node C is removed and D become the child of
#'   A. If \code{mergeSingle = FALSE}, the node C is kept. Another example: the
#'   tree has a structure as A-[B,C-[D,E]] (The root is A, and A has two
#'   children B, and C. C has two children D, and E.) If mergeSingle = TRUE,
#'   nothing changes because there is no internal node with only one child.
#' @export
#' @return A phylo object.
#'
#' @examples
#' library(ggtree)
#' data(tinyTree)
#' ggtree(tinyTree) +
#'     geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'     geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7) +
#'     geom_hilight(node = 18)
#'
#' # remove the blue branch
#' NT <- pruneTree(tree = tinyTree, rmLeaf = c(4, 5),
#'                 mergeSingle = TRUE)
#'
#' ggtree(NT) +
#'     geom_text2(aes(label = label), color = "darkorange",
#'                hjust = -0.1, vjust = -0.7)
#'

pruneTree <- function(tree, rmLeaf, mergeSingle = TRUE){

    # Use the node number
    if (is.character(rmLeaf)) {
        rmLeaf <- transNode(tree = tree, node = rmLeaf)
    }

    # edges & paths
    ed <- tree$edge
    mt <- matTree(tree = tree)

    # The path to be deleted
    mi <- match(rmLeaf, mt[, "L1"])
    nmt <- mt[-mi,]

    # find internal nodes that exist only in one path
    nodeI <- unique(ed[, 1])

    # If an internal node has only one child, merge that node. For example, a
    # tree has structure as A-[B,C-D] (The root is A, and A has two children B,
    # and C. C has a child D). If mergeSingle = TRUE, the node C is removed and
    # D become the child of A. If mergeSingle = FALSE, the node C is kept.
    # Another example: the tree has a structure as A-[B,C-[D,E]] (The root is A,
    # and A has two children B, and C. C has two children D, and E.) If
    # mergeSingle = TRUE, nothing changes because there is no internal node only
    # with one child.

    if (mergeSingle){
        locI <- lapply(nodeI, FUN = function(x){
            which(nmt == x, arr.ind = TRUE)
        })
        li <- lapply(locI, FUN = function(x) {
            nrow(x) == 1})
        li <- unlist(li)
        locII <- do.call(rbind, locI[li])

        nmt[locII] <- NA
    }

    edm <- lapply(seq_len(nrow(nmt)), FUN = function(y) {
        x <- nmt[y, ]
        xx <- x[!is.na(x)]
        rx <- rev(xx)
        lx <- length(xx)
        rxx <- cbind(rx[seq_len(lx-1)],
                     rx[setdiff(seq_len(lx), 1)])
        return(rxx)
    })

    edm <- do.call(rbind, edm)
    rownames(edm) <- NULL
    edo <- edm[!duplicated(edm), ]


    xx <- unique(sort(edo))
    mat <- cbind(old = xx, new = seq_along(xx)) # the pair (old - new )
    edn <- apply(edo, 2, FUN = function(x) {
        ind <- match(x, mat[, "old"])
        mat[ind, "new"]
    })

    # tree labels
    tipLab <- tree$tip.label
    nodeLab <- tree$node.label
    nodeA <- sort(unique(as.vector(ed)))
    names(nodeA) <- c(tipLab, nodeLab)

    tip.new <- sort(setdiff(edn[, 2], edn[, 1]))
    intNode.new <- sort(unique(edn[, 1]))

    tipLab.new <- nodeA[mat[tip.new, "old"]]
    nodeLab.new <- nodeA[mat[intNode.new, "old"]]
    rootNode <- setdiff(ed[, 1], ed[, 2])

    # if (length(tree$edge.length) < length(nodeA)) {
    #     len <- tree$edge.length
    #     ll <- cbind(len = len,
    #                 nodeNum = setdiff(nodeA, rootNode))
    #     len.new <- ll[ll[, "nodeNum"] %in% mat[, "old"], "len"]
    # } else{
    #     len <- tree$edge.length
    #     ll <- cbind(len = len, nodeNum = nodeA)
    #     len.new <- ll[ll[, "nodeNum"] %in% mat[, "old"], "len"]
    #
    # }

    len <- tree$edge.length
    len.new <- apply(edo, 1, FUN = function(x){
        distNode(tree = tree, node = x)
    })

    br <- list(edge = edn, tip.label = names(tipLab.new),
               edge.length = len.new, Nnode = length(intNode.new),
               node.label = names(nodeLab.new), nodeNumPair = mat)
    attr(br, "class") <- attr(tree, "class")
    attr(br, "order") <- attr(tree, "order")

    return(br)


}



