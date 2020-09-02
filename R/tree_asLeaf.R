#' change internal nodes to leaf nodes
#'
#' \code{asLeaf} updates a \code{phylo} tree by changing the specified internal
#' nodes to leaf nodes. In other words, the descendant nodes of the specified
#' internal nodes are removed.
#'
#' @param tree A phylo object.
#' @param node A numeric or character vector. It specifies internal nodes that
#'   are changed to leaves via their node labels or numbers.
#' @importFrom methods is
#' @export
#' @return A phylo object.
#'
#' @examples
#' library(ggtree)
#' data(tinyTree)
#' ggtree(tinyTree, ladderize = FALSE) +
#'     geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7) +
#'     geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7) +
#'     geom_hilight(node = 18) +
#'     geom_point2()
#'
#' # remove the blue branch
#' NT1 <- asLeaf(tree = tinyTree, node = 16)
#'
#' ggtree(NT1, ladderize = FALSE) +
#'     geom_text2(aes(label = label), color = "darkorange",
#'                hjust = -0.1, vjust = -0.7) +
#'     geom_point2()
#'
#' # if mergeSingle = TRUE, the node (Node_17) is removed.
#' NT2 <- asLeaf(tree = tinyTree, node = c(15, 13))
#' # or use the ape::drop.tip
#' # NT3 <-  ape::drop.tip(phy = tinyTree, tip = 4:5)
#' # all.equal(NT2, NT3)
#'
#' ggtree(NT2, ladderize = FALSE) +
#'     geom_text2(aes(label = label), color = "darkorange",
#'                hjust = -0.1, vjust = -0.7) +
#'     geom_point2()
#'
asLeaf <- function(tree, node) {

  # ===========check inputs =================
  if (!is(tree, "phylo")) {
    stop("tree should be a phylo object")
  }

  isTip <- isLeaf(tree = tree, node = node)
  if (any(isTip)) {
    stop("Leaf node(s) found:", node[isTip])
  }

  if (!is(node, "numeric") & !is(node, "character")) {
    stop("node should a numeric or character vector")
  }
  # ============remove descendants ============
  # set descendant nodes as NA
  if (is(node, "character")) {
    node <- transNode(tree = tree, node = node)
  }
  mat <- matTree(tree = tree)

  ind <- lapply(node, FUN = function(x) {
    which(mat == x, arr.ind = TRUE)
  })

  ind <- do.call(rbind, ind)

  rnd <- apply(ind, 1, FUN = function(x) {
    xx <- cbind(row = rep(x["row"], x["col"]-1),
                col = seq_len(x["col"]-1))
    # return data.frame, since matrix can get coerced
    # to vector, if length of results are all equal
    return(data.frame(xx))
  })
  rnd <- do.call(rbind, rnd)
  colnames(rnd) <- NULL
  rownames(rnd) <- NULL
  rnd <- as.matrix(rnd)
  mat[rnd] <- NA

  # move NA to the end for each row
  nat <- apply(mat, 1, FUN = function(x) {
    xx <- x[!is.na(x)]
    c(xx, rep(NA, length(x)-length(xx)))
  })
  nat <- t(nat)

  # remove duplicated rows
  natO <- nat[!duplicated(nat), ,drop=FALSE]

  # update node number
  natN <- natO
  natN[, 1] <- seq_len(nrow(natN))
  uu <- unique(sort(natN))
  od <- cbind(old = uu, new = seq_along(uu)) # the pair (old - new )
  natNN <- apply(natN, 2, FUN = function(x) {
    ind <- match(x, od[, "old"])
    od[ind, "new"]
  })

  # node pair (old - new)
  old <- as.vector(natO)
  new <- as.vector(natNN)
  pair <- cbind(old, new)
  pair <- pair[!duplicated(pair), , drop = FALSE]
  pair <- pair[rowSums(is.na(pair)) != ncol(pair), , drop = FALSE]

  if(nrow(pair) == 1L){
    stop("Selected root node of tree.")
  }
  
  # ==============Update phylo object =============================
  # new edge (edn)
  edo <- lapply(seq_len(nrow(natO)), FUN = function(y) {
    x <- natO[y, ]
    xx <- x[!is.na(x)]
    rx <- rev(xx)
    lx <- length(xx)
    rxx <- cbind(rx[seq_len(lx-1)],
                 rx[setdiff(seq_len(lx), 1)])
    return(rxx)
  })
  edo <- do.call(rbind, edo)
  rownames(edo) <- NULL
  edo <- edo[!duplicated(edo), , drop = FALSE]
  edn <- apply(edo, 2, FUN = function(x) {
    ind <- match(x, pair[, "old"])
    pair[ind, "new"]
  })

  # new edge length
  lenN <- apply(edo, 1, FUN = function(x){
    distNode(tree = tree, node = x)
  })

 # new node labels
  labN <- transNode(tree = tree, node = pair[, "old"])
  nodeA <- pair[, "new"]
  names(nodeA) <- labN
  nodeA <- sort(nodeA, decreasing = FALSE)

  labTN <- names(nodeA)[nodeA %in% natNN[, 1]]
  labIN <- names(nodeA)[!nodeA %in% natNN[, 1]]


  br <- list(edge = edn, tip.label = labTN,
             edge.length = lenN, Nnode = length(labIN),
             node.label = labIN)
  attr(br, "class") <- attr(tree, "class")
  attr(br, "order") <- attr(tree, "order")

  return(br)


}
