#' Translate a data frame to a phylo object
#'
#' \code{toTree} translates a data frame to a phylo object
#'
#' @param data A data frame or matrix.
#' @param column_order A vector that includes the column names of data to
#'   reorder columns of \code{data}. Default is NULL, the original order of
#'   \code{data} is kept.
#' @details The last column is used as the leaf nodes
#' @importFrom S4Vectors head tail
#' @importFrom dplyr arrange_all
#' @importFrom rlang .data
#' @importFrom ape reorder.phylo
#' @importFrom stats reorder
#' @return a phylo object
#' @author Ruizhu HUANG
#' @export
#' @examples
#' 
#' library(ggtree)
#' # Example 1:
#' taxTab <- data.frame(R1 = rep("A", 5),
#'                      R2 = c("B1", rep("B2", 4)),
#'                      R3 = paste0("C", 1:5))
#' # Internal nodes: their labels are prefixed with colnames of taxTab
#' #  e.g., R2:B2
#' tree <- toTree(data = taxTab)
#' # viz the tree           
#'  ggtree(tree) + 
#'  geom_text2(aes(label = label), color = "red", vjust = 1) + 
#'  geom_nodepoint()
#' 
#' # Example 2: duplicated rows in the 3rd and 4th rows
#' taxTab <- data.frame(R1 = rep("A", 5),
#'                      R2 = c("B1", rep("B2", 4)),
#'                      R3 = c("C1", "C2", "C3", "C3", "C4"))
#' # duplicated rows are removed with warnings
#' tree <- toTree(data = taxTab)
#' 
#' # Example 3: NA values in R2 column 
#' # results: the internal node with the label 'R2:'
#' taxTab <- data.frame(R1 = rep("A", 5),
#'                      R2 = c("B1", rep("B2", 2), NA, "B2"),
#'                      R3 = c("C1", "C2", "C3", NA, "C4"))
#' tree <- toTree(data = taxTab)
#' # viz the tree           
#'  ggtree(tree) + 
#'  geom_text2(aes(label = label), color = "red", vjust = 1) + 
#'  geom_nodepoint()
#' 
#' # Example 4: duplicated values in the leaf column (R4)
#' # Not allowed and give errors
#' # taxTab <- data.frame(R1 = rep("A", 5),
#' #                      R2 = c("B1", rep("B2", 3), "B3"),
#' #                      R3 = c("C1", "C2", "C3", "C3",NA),
#' #                      R4 = c("D1", "D2", "D3", NA, NA))
#' 
#' # Example 5: loops caused by missing values in B2-C4, B3-C4
#' taxTab <- data.frame(R1 = rep("A", 6),
#'                      R2 = c("B1", rep("B2", 4), "B3"),
#'                      R3 = c("C1", "C2", "C3", "C3", "C4", "C4"),
#'                      R4 = c("D1", "D2", "D3", "D3", "D4" , "D4"),
#'                      R5 = paste0("E", 1:6))
#'  # resolove loops before run to Tree
#'  # Suffix are adding to C4
#'  taxNew <- resolveLoop(taxTab)           
#'  tree <- toTree(data = taxNew)
#'  
#'  # viz the tree           
#'  ggtree(tree) + 
#'  geom_text2(aes(label = label), color = "red", vjust = 1) + 
#'  geom_nodepoint()


toTree <- function(data, column_order = NULL) {
    
    data <- data %>% data.frame() %>%
        mutate_if(is.factor, as.character) #%>%
        #replace(is.na(.data), "NA")
    # using piping this crashes if used with .data
    data <- replace(data, is.na(data), "NA")
    
    # remove duplicated rows
    isDup <- duplicated(data)
    if (any(isDup)) {
       warning(sum(isDup), " duplacted rows are removed")
        data <- data[!isDup, , drop = FALSE]
    }
    
    # check duplicated leaves
    vleaf <- data[, ncol(data)]
    if (anyDuplicated(vleaf)) {
        stop("Not allowed to have duplicated values in the leaf column: ", 
             colnames(data)[ncol(data)])
    }
    
    # the column order 
    if (!is.null(column_order)) {
        data <- data[, column_order, drop = FALSE]
    }
    column_order <- colnames(data)
    
    # The column should start from high level to low level (the root has the
    # highest level)
    ncat <- apply(data, 2, FUN = function(x) {
        length(unique(x))
    })
    if (is.unsorted(ncat)) {
        stop("The current column order: ", 
             paste(column_order, collapse = " "), ";
             Please order columns from high to low. 
             (High: the root; Low: leaves)")
    }
    
    data <- arrange_all(data.frame(data, stringsAsFactors = FALSE))
    if (ncat[[1]] !=1) {
        warning("The root is added with label 'ALL'")
        data <- cbind(root = rep("ALL", nrow(data)), data)
    }
    
    # decide leaf nodes
    vleaf <- data[, ncol(data)]
    numL <- seq_along(vleaf)
    names(numL) <- vleaf
    
    
    # decide internal nodes
    nc <- ncol(data) - 1
    
    datL1 <- lapply(seq_len(nc), FUN = function(x) {
        xx <- data[, x]
        nam <- colnames(data)[x]
        paste(nam, xx, sep = ":")
    })
    datL2 <- do.call(cbind, datL1)
    
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
        stop("The tree can't be built; loops detected in the tree.
             Try 'resolveLoop' to remove loops. ")
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
    treeList$Nnode <- length(numI)
    class(treeList) <- "phylo"
    treeList <- reorder(treeList)
    
    treeList
}


