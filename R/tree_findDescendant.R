#' Find descendants (or offsprings)
#'
#' \code{findDescendant} finds descendants of a node.
#'
#' @param node An internal node. It could be the node number or the node
#'   label.
#' @param tree A phylo object.
#' @param only.leaf A logical value, TRUE or FALSE. The default is TRUE. If
#'   default, only the leaf nodes in the descendant nodes would be returned.
#' @param self.include A logical value, TRUE or FALSE. The default is FALSE. If
#'   TRUE, the node specified in \strong{node} is included and the leaf node
#'   itself is returned as its descendant.
#' @param use.alias A logical value, TRUE or FALSE. The default is FALSE, and
#'   the node label would be used to name the output; otherwise, the alias of
#'   node label would be used to name the output. The alias of node label is
#'   created by adding a prefix \code{"alias_"} to the node number.
#' @export
#' @return A vector of nodes. The numeric value is the node number, and the
#'   vector name is the corresponding node label. If a node has no label, it
#'   would have NA as name when \code{use.alias = FALSE}, and have the alias of
#'   node label as name when \code{use.alias = TRUE}.
#' @author Ruizhu Huang
#'
#' @examples
#' data(tinyTree)
#'
#' library(ggtree)
#' ggtree(tinyTree) +
#' geom_text2(aes(label = node), color = "darkblue",
#'                hjust = -0.5, vjust = 0.7) +
#' geom_hilight(node = 17, fill = 'steelblue', alpha = 0.5) +
#' geom_text2(aes(label = label), color = "darkorange",
#'            hjust = -0.1, vjust = -0.7)
#'
#' (tips <- findDescendant(tree = tinyTree, node = c(17), only.leaf = TRUE))

findDescendant <- function(tree,
                   node,
                   only.leaf = TRUE,
                   self.include = FALSE,
                   use.alias = FALSE) {
    
    if (!inherits(tree, "phylo")) {
        stop("tree: should be a phylo object")
    }
    
    if (anyDuplicated(node)) {
        warning("duplicated values are found in input 'node'")
    }
    
    if (!(is.character(node) |
          is.numeric(node) |
          is.integer(node))) {
        stop("The argument (node) should be character or numeric")
    }
    # the edge matrix
    mat <- tree$edge
    matN <- matTree(tree = tree)
    
    if (is.character(node)) {
        numA <- convertNode(tree = tree, node = node,
                          use.alias = TRUE,
                          message = FALSE)
    } else {
        numA <- node
        isOut <- !numA %in% mat
        if (any(isOut)) {
            stop("Node ", numA,
                 " can't be found in the ",
                 deparse(substitute(tree)), "\n")
        }
        
    }
    
    loc <- lapply(seq_len(ncol(matN)), FUN = function(x){
        xx <- matN[, x]
        xm <- match(xx, numA)
        # which elements in xx could be matched to numA
        xi <- which(!is.na(xm))
        # which elements in numA
        mi <- xm[!is.na(xm)]
        
        mm <- cbind("row" = xi,
                    "col" = rep(x, length(xi)),
                    "node" = numA[mi])
        
    })
    
    moc <- do.call(rbind, loc)
    
    
    # separate leaf and internal nodes
    mocL <- moc[moc[, "col"] == 1, , drop = FALSE]
    mocI <- moc[moc[, "col"] != 1, , drop = FALSE]
    
    if (only.leaf) {
        if (!self.include) {
            mocL<- NULL
        }
        if (nrow(mocI)) {
            mocI[, "col"] <- 1
        }
        
    } else {
        if (!self.include) {
            mocL <- NULL
            if (nrow(mocI)) {
                mocI[, "col"] <- mocI[, "col"] - 1
            }
        }
        ll <- lapply(mocI[, "col"],
                     FUN = function(x) {
                         seq(from = 1, to = x , by = 1)
                     })
        mocII <- cbind("row" = rep(mocI[, "row"], mocI[, "col"]),
                       "col" = unlist(ll),
                       "node" = rep(mocI[, "node"], mocI[, "col"]))
        mocI <- mocII
    }
    
    out <- vector("list", length(numA))
    if (is.null(mocL)) {
        mocC <- mocI
    } else {
        mocC <- rbind(mocL, mocI)
    }
    
    # descendants: get, remove duplicates and sort
    desd <- cbind("found" = matN[mocC[, c("row", "col"), drop = FALSE]],
                  "node" = mocC[, "node"])
    desd <- desd[!duplicated(desd), , drop = FALSE]
    od <- order(desd[, "node"], desd[, "found"], decreasing = FALSE)
    desd <- desd[od, , drop = FALSE]
    
    # split according to the parent node
    parent <- factor(desd[, "node"],
                     levels = unique(desd[, "node"]))
    dList <- split(desd[, "found"], f = parent)
    
    # output result in the order as the input node
    #o <- match(unique(desd[, "node"]), numA)
    # out[o] <- dList
    o <- match(numA, unique(desd[, "node"]))
    out <- dList[o]
    names(out) <- convertNode(tree = tree, node = numA,
                            use.alias = use.alias)
    return(out)
}