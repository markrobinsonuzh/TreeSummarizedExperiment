#' Translate a data frame to a phylo object
#'
#' \code{toTree} translates a data frame to a phylo object
#'
#' @param data A data frame or matrix.
#' @param cache A logical value, TRUE or FALSE. The default is FALSE. If TRUE,
#'   the output `phylo` has 6 elements (edge, tip.label, edge.length, Nnode,
#'   node.label, and cache). The \strong{cache} is a list that has the length
#'   equals to the number of internal node, and each of its element stores the
#'   descendant leaves.
#'
#' @details The last column is used as the leaf nodes if the row names of the
#'   input \code{data} are not availabel; otherwise, the row names is used as
#'   the leaf nodes.
#' @importFrom utils head tail
#' @importFrom dplyr arrange_all
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
#' R2 = c("B1", rep("B2", 2), NA, "B2"),
#' R3 = c("C1", "C2", "C3", NA, "C4"))
#'
#' tree <- toTree(data = taxTab)
#'


toTree <- function(data, cache = FALSE) {

    # arrange data
    data <- as.matrix(data)
    data <- apply(data, 2, FUN = function(x) {
        x[x == "NA"] <- NA
        return(x)
    })

    data <- arrange_all(data.frame(data, stringsAsFactors = FALSE))
    # input NA with the value in the previous level
    if (any(is.na(data))) {
        cn <- colnames(data)
        data <- lapply(seq_len(ncol(data)), FUN = function(x) {
            #cat(x, "\n")
            xx <- data[, x]
            ii <- is.na(xx)
            if (sum(ii)) {
                y <- data[ii, 1: x, drop = FALSE]
                uy <- apply(y, 1, FUN = function(x) {
                    xx <- x[!is.na(x)]
                    nax <- sum(is.na(x))
                    paste(tail(xx, 1), NA, nax, sep = ":")
                })
                #xx[ii] <- paste(uy, NA, sep = ":")
                xx[ii] <- uy
            }
            return(xx)
          })
        data <- do.call(cbind, data)
        colnames(data) <- cn
    }


    # decide the leaf
    vleaf <- rownames(data)
    if (is.null(vleaf)) {
        warning("data has no rownames; Last column is used as leaf nodes. ")
        vleaf <- data[, ncol(data)]
    }

    if (anyDuplicated(vleaf)) {
        stop("Not allow to use the same label for different leaf nodes. \n")
    }

    numL <- seq_along(vleaf)
    names(numL) <- vleaf


    # decide internal nodes
    if (identical(vleaf, data[, ncol(data)])) {
        nc <- ncol(data) - 1
    } else {
        nc <- ncol(data)
    }

    datL1 <- lapply(seq_len(nc), FUN = function(x) {
        xx <- data[, x]
        nam <- colnames(data)[x]
        paste(nam, xx, sep = ":")
    })
    datL2 <- do.call(cbind, datL1)

    #nodeI <- as.vector(datL2)
    nodeI <- as.vector(t(datL2))
    nodeI <- setdiff(nodeI, vleaf)
    numI <- length(numL) + seq_along(nodeI)
    names(numI) <- nodeI

    # all nodes
    numN <- c(numL, numI)

    # create edges
    datL3 <- cbind(datL2, vleaf)
    mat1 <- apply(datL3, 2, FUN = function(x) {
        numN[x]})

    nr <- nrow(datL3)
    lx <- lapply(seq_len(nr), FUN = function(x) {
        mx <- mat1[x, ]
        mx <- mx[!is.na(mx)]
        cx <- cbind(head(mx, -1), tail(mx, -1))
        cx})

    mat2 <- do.call(rbind, lx)
    mat3 <- mat2[!duplicated(mat2), ]


    # sort node number
    numL <- sort(numL)
    numI <- sort(numI)

    # tree
    treeList <- vector("list", 5)
    names(treeList) <- c("edge", "tip.label", "node.label",
                         "edge.length", "Nnode")
    treeList$edge <- matrix(mat3, ncol = 2)
    treeList$tip.label <- names(numL)
    treeList$node.label <- names(numI)
    treeList$edge.length <- rep(0.1, nrow(mat3))
    #treeList$edge.length <- runif(nrow(mat3))
    treeList$Nnode <- length(numI)
    class(treeList) <- "phylo"
    attr(treeList, "order") <- "cladewise"

    # keep cache
    desA2 <- findOS(tree = treeList, ancestor = numI, only.Tip = TRUE)
    desA1 <- split(numL, names(numL))
    desA <- c(desA1, desA2)

    if (cache) {treeList$cache <- desA}

    return(treeList)
}
