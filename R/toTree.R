#' Translate a data frame to a phylo object
#'
#' \code{toTree} translates a data frame to a phylo object
#'
#' @param data A data frame or matrix.
#'
#' @importFrom utils head tail
#' @return a phylo object
#' @author Ruizhu HUANG
#' @export
#' @examples
#'
#' taxTab <- data.frame(R1 = rep("A", 5),
#' R2 = c("B1", rep("B2", 4)),
#' R3 = c("C1", "C2", "C3", "C3", "C4"))
#'
#' taxTab <- data.frame(R1 = rep("A", 5),
#' R2 = c("B1", rep("B2", 4)),
#' R3 = c(NA, "C2", "C3", NA, "C4"))
#'
#' tree <- toTree(data = taxTab)
#'


toTree <- function(data, cache = FALSE) {

    # input NA with the value in the previous level
    if (any(is.na(data))) {
        data <- apply(data, 1, FUN = function(x) {
            xx <- x[!is.na(x)]
            x[is.na(x)] <- paste(tail(xx, 1), NA, sep = " - ")
            return(x)
        })
        data <- t(data)
    }


    # add level
    datL1 <- lapply(seq_len(ncol(data)), FUN = function(x) {
        xx <- data[, x]
        nam <- colnames(data)[x]
        paste(nam, "-", xx)
    })
    datL2 <- do.call(cbind, datL1)

    # leaf nodes
    vleaf <- rownames(data)

    if (is.null(vleaf)) {
        vleaf <- datL2[, ncol(datL2)]
    }
    if (is.null(vleaf)) {
        stop("Not allow to use the same label for different leaf nodes. \n")
    }

    if (any(vleaf != datL2[, ncol(datL2)])) {
        datL2 <- cbind(datL2, vleaf)
    }

    numL <- seq_along(vleaf)
    names(numL) <- vleaf

    # internal nodes
    nr <- nrow(datL2)
    lx1 <- lapply(seq_len(nr), FUN = function(x) {
        tx <- unlist(datL2[x, ])
        tx <- tx[!is.na(tx)]
        setdiff(tx, vleaf)
    })
    vlx1 <- unlist(lx1)
    vlx1 <- vlx1[!duplicated(vlx1)]
    numI <- length(numL) + seq_along(vlx1)
    names(numI) <- vlx1

    # all nodes
    numN <- c(numL, numI)

    mat1 <- apply(datL2, 2, FUN = function(x) {
        numN[x]})

    lx2 <- lapply(seq_len(nr), FUN = function(x) {
        mx <- mat1[x, ]
        mx <- mx[!is.na(mx)]
        cx <- cbind(head(mx, -1), tail(mx, -1))
        cx})

    mat2 <- do.call(rbind, lx2)
    mat3 <- mat2[!duplicated(mat2), ]


    # sort node number
    numL <- sort(numL)
    numI <- sort(numI)

    # tree
    treeList <- vector("list", 5)
    names(treeList) <- c("edge", "tip.label", "node.label",
                         "edge.length", "Nnode")
    treeList$edge <- mat3
    treeList$tip.label <- names(numL)
    treeList$node.label <- names(numI)
    treeList$edge.length <- rep(0.1, nrow(mat3))
    treeList$Nnode <- length(numI)
    class(treeList) <- "phylo"

    # keep cache
    node <- unique(as.vector(datL2))
    desA <- lapply(node, FUN = function(x) {
        xi <- which(datL2== x, arr.ind = TRUE)
        ri <- xi[, "row"]
        lvs <- datL2[ri, ncol(datL2)]
        ln <- transNode(tree = treeList, input = lvs)
        names(ln) <- lvs
    })
    names(desA) <- node

    treeList$cache <- desA

    return(treeList)
}
