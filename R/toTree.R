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
#'
#' tree <- toTree(data = taxTab)
#'


toTree <- function(data) {

    # input NA with the value in the previous level
    data <- apply(data, 1, FUN = function(x) {
        xx <- x[!is.na(x)]
        x[is.na(x)] <- tail(xx, 1)
        return(x)
    })
    data <- t(data)

    # add level
    datL1 <- lapply(seq_len(ncol(data)), FUN = function(x) {
        xx <- data[, x]
        nam <- colnames(data)[x]
        paste(nam, "-", xx)
    })
    datL2 <- do.call(cbind, datL1)

    # leaf nodes
    # leaf <- apply(datL2, 1, FUN = function(x) {
    #     #tail(x[x != "Unclassified"], 1)
    #     tail(x[! is.na(x)], 1)
    # })
    # vleaf <- as.vector(leaf)
    # vleaf <- unique(vleaf)
    vleaf <- rownames(data)
    if (any(vleaf != datL2[, ncol(datL2)])) {
        datL2 <- cbind(datL2, vleaf)
    }

    # leaves
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
    return(treeList)
}
