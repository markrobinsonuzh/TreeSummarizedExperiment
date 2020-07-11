#' list leaf nodes that are the descendants of at least one specified node
#'
#' \code{unionLeaf} list the leaf nodes that are the desendants of (at least one)
#' specified nodes.
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
#' u1 <- unionLeaf(tree = tinyTree, node = c(19, 17))
#' u2 <- unionLeaf(tree = tinyTree, node = c(19, 17, 7))
#' (u3 <- unionLeaf(tree = tinyTree, node = c(11, 17, 7)))
unionLeaf <- function(tree, node) {

  # ===========check inputs =================
  if (!is(tree, "phylo")) {
    stop("tree should be a phylo object")
  }

  if (!is(node, "numeric") & !is(node, "character")) {
    stop("node should a numeric or character vector")
  }
  # ============remove descendants ============
  if (is(node, "character")) {
    node <- transNode(tree = tree, node = node)
  }

  # all paths that connecting leaves and the root
  mat <- matTree(tree = tree)

  loc <- lapply(seq_len(ncol(mat)), FUN = function(x){
    xx <- mat[, x]
    xm <- match(xx, node)

    # which elements in xx could be matched to numA
    xi <- which(!is.na(xm))
    return(xi)
    })
  loc <- unlist(loc)
  loc <- unique(loc)

  moc <- cbind("row" = loc, "col" = rep(1, length(loc)))

  out <- sort(mat[moc])
  return(out)
  }
