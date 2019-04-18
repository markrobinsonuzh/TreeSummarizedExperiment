
#' Perform data aggregations based on the available tree structures
#'
#' \code{aggValue} aggregates values on the leaf nodes of a tree to a specific
#' arbitrary level of the tree. The level is specified via the nodes of the
#' tree. Users could decide on which dimension (row or column) and how should
#' the aggregation be performed.
#'
#' @param x A \code{TreeSummarizedExperiment} object or a matrix. If the latter
#'   is given, the tree structure is required (more details in \code{rowTree}
#'   and \code{colTree}).
#' @param rowLevel A numeric or character vector. The default is NULL. It
#'   provides the level on the tree that the aggregation is performed to. The
#'   aggregation is on the row dimension of matrix (or the \code{assay} tables).
#' @param colLevel A numeric or character vector. The default is NULL. It
#'   provides the level on the tree that the aggregation is performed to. The
#'   aggregation is on the column dimension of matrix (or the \code{assay}
#'   tables).
#' @param FUN A function to be applied on the aggregation. It's similar to the
#'   \code{FUN} in \code{\link[base]{apply}}
#' @param message A logical value. The default is TRUE. If TRUE, it will print
#'   out the running process.
#' @import SingleCellExperiment
#' @importFrom utils flush.console
#' @return A \code{TreeSummarizedExperiment} object or a
#'   \code{matrix}. The output has the same class of the input \code{x}.
#' @export
#'
#' @include allClass.R
#' @author Ruizhu HUANG
#' @seealso \code{\link{TreeSummarizedExperiment}}
#' @examples
#'
#' # assays data
#' set.seed(1)
#' toyTable <- matrix(rnbinom(20, size = 1, mu = 10), nrow = 5)
#' colnames(toyTable) <- paste(rep(LETTERS[1:2], each = 2),
#'                             rep(1:2, 2), sep = "_")
#' rownames(toyTable) <- paste("entity", seq_len(5), sep = "")
#'
#' toyTable
#'
#' # the column data
#' colInf <- DataFrame(gg = c(1, 2, 3, 3),
#'                     group = rep(LETTERS[1:2], each = 2),
#'                     row.names = colnames(toyTable))
#' colInf
#'
#' # the toy tree
#' library(ape)
#' set.seed(4)
#' treeC <- rtree(4)
#' treeC$node.label <- c("All", "GroupA", "GroupB")
#'
#' library(ggtree)
#' ggtree(treeC, size = 2) +
#'    geom_text2(aes(label = node), color = "darkblue",
#'            hjust = -0.5, vjust = 0.7, size = 6) +
#'     geom_text2(aes(label = label), color = "darkorange",
#'                hjust = -0.1, vjust = -0.7, size = 6)
#'
#' colInf$nodeLab <- treeC$tip.label
#' tse <- TreeSummarizedExperiment(assays = list(toyTable),
#'                                 colData = colInf,
#'                                 colTree = treeC)
#'
#' aggCol <- aggValue(x = tse, colLevel = c("GroupA", "GroupB"),
#' FUN = sum)
#'
#' assays(aggRow)[[1]]
aggValue <- function(x, rowLevel = NULL, colLevel = NULL,
                    FUN = sum, message = FALSE) {
    ## --------------------- check input before run ---------------------------
    ## If x is a TreeSummarizedExperiment, rowTree and colTree should be set to
    ## NULL
    if (!is(x, "TreeSummarizedExperiment")) {
        stop("x should be a TreeSummarizedExperiment")

    }

    ## The provided rowLevel or colLevel should be a numeric or character vector
    isR <- (is.character(rowLevel) | is.numeric(rowLevel) |
             is.null(rowLevel))
    isC <- (is.character(colLevel) | is.numeric(colLevel) |
            is.null(colLevel))
    if (!isR) {
        stop("rowLevel should be a numeric or character vector. \n")
    }
    if (!isC) {
        stop("colLevel should be a numeric or character vector. \n")
    }

    ## ---------------------- get all data ready ------------------------------
    if (message) {message("Preparing data... ")}

    ## The assay tables
    mat <- assays(x)

    ## The trees
    rTree <- rowTree(x)
    cTree <- colTree(x)

    ## Indicate whether aggregation should be on rows or columns
    if (is.null(rowLevel)) {
        onRow <- FALSE
    } else {
        if (is.null(rTree)) {
            stop("The tree on rows doesn't exist.")
        }

        onRow <- TRUE
        if (message) {
            message("Perform aggregation on the row dimension... ")
        }
       }

    if (is.null(colLevel)) {
        onCol <- FALSE
    } else {
        if (is.null(cTree)) {
            stop("The tree on columns doesn't exist.")
        }

        onCol <- TRUE
        if (message) {
            message("Perform aggregation on the column dimension... ")
        }
    }

    ## -------------------- aggregation on row dimension ----------------------
    if (onRow) {
        rD <- rowData(x)
        rL <- rowLink(x)
        outR <- .aggFun(tree = rTree,
                        assayTab = mat,
                        dimData = rD,
                        linkData = rL,
                        level = rowLevel,
                        FUN = FUN,
                        message = message)
        nrD <- outR$newDD
        mat <- outR$dataTab
    }else{
        nrD <- rowData(x)
        nrD$nodeLab <- rowLink(x)$nodeLab
    }

    ## -------------------- aggregation on column dimension ----------------
    if (onCol) {
        cD <- colData(x)
        cL <- colLink(x)
        outC <- .aggFun(tree = cTree,
                        assayTab = lapply(mat, t),
                        dimData = cD,
                        linkData = cL,
                        level = colLevel,
                        FUN = FUN)
        ncD <- outC$newDD
        mat <- lapply(outC$dataTab, t)
    } else {
        ncD <- colData(x)
        ncD$nodeLab <- colLink(x)$nodeLab
    }

    # create the new TreeSummarizedExperiment object
    out <- TreeSummarizedExperiment(assays = mat,
                                    rowTree = rTree,
                                    colTree = cTree,
                                    rowData = nrD,
                                    colData = ncD)

    return(out)
}


.aggFun <- function(tree, assayTab, dimData, linkData,
                    level = NULL, FUN, message = FALSE) {

    # nodeNum
    numR <- linkData$nodeNum

    # All node numbers on the tree
    ed <- tree$edge
    numAR <- unique(as.vector(ed))

    # The aggregation level
    if (is.null(level)) {level <- numR}
    if (is.character(level)) {
        level <- transNode(tree = tree, node = level,
                           use.alias = FALSE, message = FALSE)
    }

    # The descendants of the aggregation level
    if (is.null(tree$cache)) {
        desR <- findOS(tree = tree, node = level,
                       only.leaf = TRUE, self.include = TRUE)
        names(desR) <- transNode(tree = tree, node = level,
                                 use.alias = TRUE, message = FALSE)
    } else {
        desR <- tree$cache[level]
        names(desR) <- transNode(tree = tree, node = level,
                                 use.alias = TRUE, message = FALSE)
    }

    # Find the rows of the descendants
    idR <- lapply(desR, FUN = function(x) {
        numR %in% x
    })

    miR <- lapply(desR, FUN = function(x){
        x[!x %in% numR]
    })
    miR <- unique(unlist(miR))
    if (length(miR)) {
        warning(length(miR), "leaves couldn't be found from the assay table.
                You might want to clean the tree before aggregation \n")
    }

    ## Perform the aggregation
    # on the assay tables

    if (message) {
        message("Working on the assays table... ")}

    outR <- vector("list", length(assayTab))
    names(outR) <- names(assayTab)

    ll <- numeric()
    for (i in seq_along(assayTab)) {
        mi <- assayTab[[i]]
        oi <- lapply(seq_along(idR), FUN = function(x) {

            if (message) {
                message(x, " out of ", length(idR),
                        " finished", "\r", appendLF = FALSE)
                flush.console()
            }

            mx <- mi[idR[[x]], , drop = FALSE]
            ax <- apply(mx, 2, FUN = FUN)
            rx <- rbind(ax)
            rownames(rx) <- rep(names(idR)[x], nrow(rx))
            return(rx)
        })

        roi <- do.call(rbind, oi)
        outR[[i]] <- roi

        # record the output length
        if (i == 1) {
            lo <-  lapply(oi, nrow)
            ll <- unlist(lo)
        }
    }

    # on the dimData (rowData/colData)
    if (message) {
        message("Working on the row/column data... ")}


    nc <- ncol(dimData)
    newDD <- dimData[rep(1, sum(ll)), ]
    for (i in seq_len(nc)) {
        ri <- lapply(seq_along(idR), FUN = function(x) {
            if (message) {
                message(x, " out of ", length(idR),
                        " finished", "\r", appendLF = FALSE)
                flush.console()
            }

            xx <- idR[[x]]
            rx <- dimData[xx, i]
            ru <- unique(rx)
            if (length(ru) > 1) {
                ui <- NA
            } else {
                ui <- ru
            }

            ul <- rep(ui, ll[x])
            return(ul)
        })
        newDD[, i] <- unlist(ri)

    }
    newDD$nodeLab <- rep(names(idR), ll)

    # use the same row names as outR
    rownames(newDD) <- rep(names(idR), ll)

    out <- list(dataTab = outR, newDD = newDD)
    return(out)
}









