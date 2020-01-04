#' Translate a data frame to a phylo object
#'
#' \code{toTree} translates a data frame to a phylo object
#'
#' @param data A data frame or matrix.
#' @param edges A logical value, TRUE or FALSE. The default is FALSE. If TRUE,
#'   the input data should be a matrix of edges (A matrix with two columns). The
#'   first column is the parent and the second is the child. If FALSE,
#'   \code{data} should be a taxonomic table (see examples).
#'
#' @details The last column is used as the leaf nodes
#' @importFrom S4Vectors head tail
#' @importFrom dplyr arrange_all
#' @importFrom ape reorder.phylo
#' @importFrom stats reorder
#' @return a phylo object
#' @author Ruizhu HUANG
#' @export
#' @examples
#' # data is a taxonomic table
#' taxTab <- data.frame(R1 = rep("A", 5),
#' R2 = c("B1", rep("B2", 4)),
#' R3 = c("C1", "C2", "C3", "C3", "C4"))
#'
#' taxTab <- data.frame(R1 = rep("A", 5),
#' R2 = c("B1", rep("B2", 2), NA, "B2"),
#' R3 = c("C1", "C2", "C3", NA, "C4"))
#'
#' tree <- toTree(data = taxTab, edges = FALSE)
#' 
#' 
#' # data is a matrix of edges
#' # the matrix should be in character
#'  mat <- tree$edge
#'  mat <- apply(mat, 2, function(x){
#'           transNode(tree, node = x)})
#'  tree2 <- toTree(mat, edges = TRUE)
#'  

toTree <- function(data, edges = FALSE) {
    if (edges) {
        out <- .toTree2(data)
    } else {
        out <- .toTree1(data)
    }
}

.toTree2 <- function(data) {
    data <- as.matrix(data)
    if (!is.character(as.vector(data))) {
        stop("data should be a character matrix")
    }
    
    # leaf nodes
    leaf <- setdiff(data[, 2], data[, 1])
    
    # internal nodes
    inNode <- setdiff(data[, 1], leaf)
    
    # all nodes
    nodeLab <- c(leaf, inNode)
    node <- seq_along(nodeLab)
    names(node) <- nodeLab
    
    # edges
    edges <- apply(data, 2, FUN = function(x) {node[x]})
    
    # tree
    tree <- list(
        edge = edges,
        tip.label= leaf,
        Nnode = length(node),
        node.label = inNode)
    
    class(tree)<-"phylo"
    tree <- reorder.phylo(tree)
    
    return(tree)
}
.toTree1 <- function(data) {

    # ========== check ========================
    # arrange data
    data <- as.matrix(data)
    data <- apply(data, 2, FUN = function(x) {
        x[x == "NA"] <- NA
        return(x)
    })

    # The column should start from high level to low level (the root has the
    # highest level)
    ncat <- apply(data, 2, FUN = function(x) {
        length(unique(x))
    })
    if (is.unsorted(ncat)) {
        stop("Columns should be ordered from the high to low (leaf level)")
    }
    
    data <- arrange_all(data.frame(data, stringsAsFactors = FALSE))
    if (ncat[[1]] !=1) {
        warning("The root is added with label 'ALL'")
        data <- cbind(root = rep("ALL", nrow(data)), data)
    }
    # input NA with the value in the previous level
    if (any(is.na(data))) {
        cn <- colnames(data)
        data <- lapply(seq_len(ncol(data)), FUN = function(x) {
            #cat(x, "\n")
            xx <- data[, x]
            ii <- is.na(xx)
            if (sum(ii)) {
                y <- data[ii, seq_len(x), drop = FALSE]
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
    #vleaf <- rownames(data)
    # if (is.null(vleaf)) {
    #     warning("data has no rownames; Last column is used as leaf nodes. ")
    #     vleaf <- data[, ncol(data)]
    # }

    vleaf <- data[, ncol(data)]
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

    tt <- table(mat3[, 2])
    if (any(tt > 1)) {
        bt <- tt[tt>1]
        st <- as.numeric(names(bt))
        lt <- gsub(pattern = ".*:",  replacement = "",
                   names(numN[st]))
        ct <- gsub(pattern = ":.*",  replacement = "",
                    names(numN[st]))
        dt <- data.frame(node = st,
                         Freq = as.vector(bt),
                         column = ct,
                         value = lt)
        stop("The tree can't be built; loops detected in the tree.")
    }

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
    treeList <- reorder(treeList)

    return(treeList)
}
