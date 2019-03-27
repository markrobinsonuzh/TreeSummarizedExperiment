#' Transform a phylo object into a matrix.
#'
#' \code{matTree} transforms a phylo tree into a matrix. The entry of the matrix
#' is node number. Each row represents a path connecting a leaf node and the
#' root. The columns are arranged in the order as the path passing the nodes to
#' reach the root.
#'
#' @param tree A phylo object
#' @export
#' @return A matrix
#' @author Ruizhu Huang
#'
#' @examples
#' library(ggtree)
#'
#' data(tinyTree)
#' ggtree(tinyTree, branch.length = 'none') +
#'  geom_text2(aes(label = node))
#'
#'
#' # each row of the matrix representing a path.
#' # the first column is leaf nodes; the last non-NA value in a row is the root
#' mat <- matTree(tree = tinyTree)
#'
matTree <- function(tree) {

    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }


    # the edge matrix
    mat <- tree$edge

    # the leaves
    L1 <- setdiff(mat[, 2], mat[, 1])

    # each path connects a tip with the root.
    # each path is stored as a row of matN
    # the first column stores the tips (or leaves)
    matN <- cbind(L1)
    i <- 1
    repeat {
        li <- mat[match(matN[, i], mat[, 2]), 1]
        ll <- length(unique(li[!is.na(li)]))
        if (ll == 0) {
            break
        }
        matN <- cbind(matN, li)
        i <- i + 1
    }
    rownames(matN) <- NULL
    colnames(matN) <- paste("L", seq_len(ncol(matN)), sep = "")

    return(matN)
}

