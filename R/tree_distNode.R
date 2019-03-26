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
#' # remove the blue branch
#' distNode(tree = tinyTree, node = c(10, 11))
#'



distNode <- function(tree, node) {

    if (is.character(node)) {
        node <- transNode(tree = tree, node = node)
    }

    sNode <- shareNode(tree = tree, node = node)

    ed <- tree$edge
    mt <- matTree(tree = tree)
    wi <- apply(mt, 1, FUN = function(x) {
        nn <- c(node, sNode)
        nn <- unique(nn)
        sum(nn %in% x) == 2
    })


    # edges connect two nodes
    if (sum(wi) > 2) {
        wi <- which(wi)[1]
        smt1 <- mt[wi, ]
        ss <- smt1[!is.na(smt1)]
        rs <- rev(ss)
        smt2 <- rs[rs %in% node]
        smt2 <- rbind(smt2)

    } else {
        smt1 <- mt[wi, , drop = FALSE]
        smt2 <- lapply(seq_len(nrow(smt1)), FUN = function(x){
            xx <- smt1[x, ]
            wx <- which(xx == sNode)
            x1 <- xx[seq_len(wx)]
            x2 <- rev(x1)
            rr <- seq_len(length(x2)-1)
            cbind(x2[rr], x2[rr+1])
        })
        smt2 <- do.call(rbind, smt2)
    }

    fi <- match(data.frame(t(smt2)), data.frame(t(ed)))

    len <- tree$edge.length

    dd <- sum(len[fi])

    return(dd)
}
