#' Convert a data frame to a phylo object
#'
#' \code{asPhylo} converts a data frame to a phylo object. Compared to
#' \code{toTree}, \code{asPhylo} allows the output tree to have different number
#' of nodes in paths connecting leaves to the root.
#'
#' @param data A data frame or matrix.
#' @param column_order A vector that includes the column names of data to
#'   reorder columns of \code{data}. Default is NULL, the original order of
#'   \code{data} is kept.
#' @param asNA This specifies strings that are considered as NA
#' @details The last column is used as the leaf nodes
#' @importFrom treeio as.phylo
#' @importFrom rlang .data
#' @return a phylo object
#' @author Ruizhu HUANG
#' @export
#' @examples
#' 
#' library(ggtree)
#' 
#' # Example 0:
#' taxTab <- data.frame(R1 = rep("A", 5),
#'                      R2 = c("B1", rep("B2", 4)),
#'                      R3 = paste0("C", 1:5))
#' # Internal nodes: their labels are prefixed with colnames of taxTab
#' #  e.g., R2:B2
#' taxTree <- asPhylo(data = taxTab)
#' ggtree(taxTree) + 
#'  geom_text2(aes(label = label), color = "red", vjust = 1) + 
#'  geom_nodepoint()
#'  
#'  # (Below gives the same output as toTree)
#'  taxTab$R1 <- paste0("R1:", taxTab$R1)
#'  taxTab$R2 <- paste0("R2:", taxTab$R2)
#'  taxTree <- asPhylo(data = taxTab)
#' # viz the tree           
#'  ggtree(taxTree) + 
#'  geom_text2(aes(label = label), color = "red", vjust = 1) + 
#'  geom_nodepoint()
#'  
#' # Example 1
#' df1 <- rbind.data.frame(c("root", "A1", "A2", NA),
#'                         c("root", "B1", NA, NA))
#' colnames(df1) <- paste0("L", 1:4)
#' tree1 <- asPhylo(df1)
#' 
#' ggtree(tree1, color = "grey") +
#' geom_nodepoint() +
#'    geom_text2(aes(label = label), angle = 90,
#'                color = "red", vjust = 2,
#'                size = 4)
#'                
#' # Example 2
#' df2 <- data.frame(Group_1 = rep("Root", 11),
#'                   Group_2 = rep(c(13, 21), c(9, 2)),
#'                   Group_3 = rep(c(14, 18, "unknown"), c(5, 4, 2)),
#'                   Group_4 = rep(c(15, "unknown", 19, "unknown"), c(4, 1, 3, 3)),
#'                   Group_5 = rep(c(16, "unknown", 20, "unknown"), c(3, 2, 2, 4)),
#'                   Group_6 = rep(c(17, "unknown"), c(2, 9)),
#'                   LEAF = 1:11)
#'                   
#' tree2 <- asPhylo(df2, asNA = "unknown")
#' 
#' ggtree(tree2, color = "grey") +
#' geom_nodepoint() +
#'    geom_text2(aes(label = label), angle = 90,
#'                color = "red", vjust = 2,
#'                size = 4) 
#'                
#' # Example 3
#' df3 <- df2
#' df3[10:11, 3] <- ""                                            
#'  
#' tree3 <- asPhylo(df3, asNA = c("unknown", ""))
#' 
#' ggtree(tree3, color = "grey") +
#' geom_nodepoint() +
#'    geom_text2(aes(label = label), angle = 90,
#'                color = "red", vjust = 2,
#'                size = 4) 

asPhylo <- function(data, column_order = NULL, asNA = NULL) {
    
    # the column order 
    if (!is.null(column_order)) {
        data <- data[, column_order, drop = FALSE]
    }
    
    # convert specified strings as NA
    if (!is.null(asNA)) {
        ind <- lapply(asNA, function(x) {
            which(data == x, arr.ind = TRUE)
        })
        ind <- do.call(rbind, ind)
        data[ind] <- NA
    }
    
    # assign nodeID
    leaf <- apply(data, 1, .capture_last_noNA)
    id <- sort(.assign_node(data = data))
    data <- .replace_with_nodeID(data = data)
    
    # move NA to the end of each row
    data <- t(apply(data, 1, .move_na_to_end))
    
    # convert to edges
    ed <- .convert_to_edges(data = data)
    ed <- data.frame(ed)
    # convert it into phylo
    out <- as.phylo(ed)
    out$tip.label <- names(id)[seq_along(leaf)]
    out$node.label <- names(id)[-seq_along(leaf)]
    return(out)
}

.move_na_to_end <- function(x) {
    na <- is.na(x)
    c(x[!na], x[na])
}

.capture_last_noNA <- function(x) {
    na <- is.na(x)
    tail(x[!na], 1)
}

.assign_node <- function(data){
    leaf <- apply(data, 1, .capture_last_noNA)
    
    # internal nodes
    vv <- as.vector(t(data))
    vv <- vv[!is.na(vv)]
    inNode <- unique(setdiff(vv, leaf))
    
    lab <- c(unique(leaf), inNode)
    nodeID <- setNames(seq_along(lab), lab)
    
    return(nodeID)
    
}

.replace_with_nodeID <- function(data){
    id <- .assign_node(data)
    apply(data, 2, function(x) {
        id[x]
    })
}

.convert_to_edges <- function(data){
    data <- t(apply(data, 1, .move_na_to_end))
    n <- seq_len(ncol(data)-1)
    ed <- lapply(n, function(i) {
        xi <- data[, i:(i+1), drop = FALSE]
        na_i <- which(is.na(xi), arr.ind = TRUE)[, "row"]
        if (length(na_i)) {
            xi <- xi[-na_i, ]
        }
        
        unname(unique(xi))
        
    })
    
    ed <- do.call(rbind, ed)
    unique(ed)
}