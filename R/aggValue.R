# Aggregation with input x as a TreeSummarizedExperiment
.aggMatrix <- function(x, rowTree = NULL, colTree = NULL,
                       rowLevel = NULL, colLevel = NULL,
                       FUN = sum, message = FALSE) {

    ## Indicate whether aggregation should be on rows or columns
    if (is.null(rowTree)) {
        if (!is.null(rowLevel)) {
            stop("The tree on rows doesn't exist.")
        }
        onRow <- FALSE
    } else {
        onRow <- TRUE
    }

    if (is.null(colTree)) {
        if (!is.null(colLevel)) {
            stop("The tree on columns doesn't exist.")
        }
        onCol <- FALSE
    } else {
        onCol <- TRUE
    }

    if (onRow) {
        mat <- list(x)
        datM <- NULL
        datL <- DataFrame(nodeLab = rownames(x),
                          nodeNum = transNode(tree = rowTree,
                                              node = rownames(x),
                                              use.alias = FALSE,
                                              message = FALSE),
                          isLeaf = isLeaf(tree = rowTree,
                                          node = rownames(x)))
        outR <- .aggFun(tree = rowTree,
                        assayTab = mat,
                        mCol = datM,
                        linkData = datL,
                        level = rowLevel,
                        FUN = FUN)
        x <- (outR$dataTab)[[1]]

        # rownames
        newRD <- as.data.frame(outR$LinkDF)
        lab <- newRD$nodeLab_alias
        if (is.null(lab)) {
            lab <- newRD$nodeLab
        }
        rownames(x) <- lab
    }

    if (onCol) {
        mat <- list(t(x))
        datM <- NULL
        datL <- DataFrame(nodeLab = colnames(x),
                          nodeNum = transNode(tree = colTree,
                                              node = colnames(x),
                                              use.alias = FALSE,
                                              message = FALSE),
                          isLeaf = isLeaf(tree = colTree,
                                          node = rownames(x)))
        outR <- .aggFun(tree = colTree,
                        assayTab = mat,
                        mCol = datM,
                        linkData = datL,
                        level = colLevel,
                        FUN = FUN)
        xx <- (outR$dataTab)[[1]]
        x <- t(xx)

        # rownames
        newRD <- as.data.frame(outR$LinkDF)
        lab <- newRD$nodeLab_alias
        if (is.null(lab)) {
            lab <- newRD$nodeLab
        }
        colnames(x) <- lab
    }

    return(x)
}

# Aggregation with input x as a matrix
.aggTSE <- function(x, rowTree = NULL, colTree = NULL,
                    rowLevel = NULL, colLevel = NULL,
                    FUN = sum, message = FALSE) {
    ## --------------------- check input before run ---------------------------
    ## If x is a TreeSummarizedExperiment, rowTree and colTree should be set to
    ## NULL
    if (is(x, "TreeSummarizedExperiment")) {
        if (!is.null(rowTree) | !is.null(colTree)) {
            stop("The input x is a TreeSummarizedExperiment. \n",
                 "rowTree and colTree should both be NULL.\n ")
        }
    }



    ## The provided rowLevel or colLevel should be a numeric or character vector
    isR <- (is.character(rowLevel) | is.numeric(rowLevel) |
                is.integer(rowLevel) | is.null(rowLevel))
    isC <- (is.character(colLevel) | is.numeric(colLevel) |
                is.integer(colLevel) | is.null(colLevel))
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
    rowTree <- treeData(x, onRow = TRUE)
    colTree <- treeData(x, onRow = FALSE)

    ## Indicate whether aggregation should be on rows or columns
    if (is.null(rowTree)) {
        if (!is.null(rowLevel)) {
            stop("The tree on rows doesn't exist.")
        }
        onRow <- FALSE
    } else {
        if (!is.null(rowLevel)) {
            onRow <- TRUE
            if (message) {
                message("Perform aggregation on the row dimension... ")
                }
        } else {
            onRow <- FALSE
        }

    }

    if (is.null(colTree)) {
        if (!is.null(colLevel)) {
            stop("The tree on columns doesn't exist.")
        }
        onCol <- FALSE
    } else {
        if (!is.null(colLevel)) {
            onCol <- TRUE
            if (message) {
                message("Perform aggregation on the column dimension... ")
            }
        } else {
            onCol <- FALSE
        }

    }



    ## -------------------- aggregation on row dimension ----------------------
    if (onRow) {
        mat <- assays(x)
        rowM <- metaCol(x, onRow = TRUE)
        rowL <- linkData(x, onRow = TRUE)
        outR <- .aggFun(tree = rowTree,
                        assayTab = mat,
                        mCol = rowM,
                        linkData = rowL,
                        level = rowLevel,
                        FUN = FUN,
                        message = message)
        newRD <- outR$LinkDF
        mat <- outR$dataTab

        # update the row names
        rn <- as.data.frame(newRD)$nodeLab
        if (any(is.na(rn))) {
            rownames(newRD) <- NULL
            mat <- lapply(mat, FUN = function(x) {
                rownames(x) <- NULL
                return(x)
            })
        } else {
            rownames(newRD) <- rn
            mat <- lapply(mat, FUN = function(x) {
                rownames(x) <- rn
                return(x)
            })
        }




    }else{
        newRD <- rowData(x)
    }

    ## -------------------- aggregation on column dimension ----------------------
    if (onCol) {
        colM <- metaCol(x, onRow = FALSE)
        colL <- linkData(x, onRow = FALSE)
        outC <- .aggFun(tree = colTree,
                        assayTab = lapply(mat, t),
                        mCol = colM,
                        linkData = colL,
                        level = colLevel,
                        FUN = FUN)
        newCD <- outC$LinkDF
        mat <- lapply(outC$dataTab, t)

        # update the row names
        rn <- as.data.frame(newCD)$nodeLab
        if (any(is.na(rn))) {
            rownames(newCD) <- NULL
            mat <- lapply(mat, FUN = function(x) {
                colnames(x) <- NULL
                return(x)
            })
        } else {
            rownames(newCD) <- rn
            mat <- lapply(mat, FUN = function(x) {
                colnames(x) <- rn
                return(x)
            })
        }


    } else {
        newCD <- colData(x)
    }

    # create the new TreeSummarizedExperiment object
    out <- TreeSummarizedExperiment(rowTree = rowTree, colTree = colTree,
                                    assays = mat, rowData = newRD,
                                    colData = newCD)

    return(out)
}


.aggFun <- function(tree, assayTab, mCol, linkData,
                    level, FUN, message = FALSE) {

    # The link information via nodeNum
    numR <- linkData$nodeNum

    # All node numbers on the tree
    ed <- tree$edge
    numAR <- unique(as.vector(ed))

    # The aggregation level
    if (is.null(level)) {level <- numAR}
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
        warning(length(miR), "leaves couldn't be found from the data table.
                The value 0 is used for the missing leaves. \n")
    }

    ## Perform the aggregation
    # on the assay tables

    if (message) {
        message("Working on the assays table... ")}

    outR <- vector("list", length(assayTab))
    names(outR) <- names(assayTab)

    ll <- numeric()
    for (i in seq_along(outR)) {
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
        # colnames(roi) <- colnames(mi)
        # rownames(roi) <- NULL
        outR[[i]] <- roi

        # record the output length
        if (i == 1) {
            lo <-  lapply(oi, nrow)
            ll <- unlist(lo)
        }
    }

    # on the row data
    if (message) {
        message("Working on the row/column data... ")}

    if (!is.null(mCol)) {
        nc <- ncol(mCol)
        rnCol <- mCol[rep(1, sum(ll)), ]
        for (i in seq_len(nc)) {
            ri <- lapply(seq_along(idR), FUN = function(x) {
                if (message) {
                    message(x, " out of ", length(idR),
                            " finished", "\r", appendLF = FALSE)
                    flush.console()
                }

                xx <- idR[[x]]
                rx <- mCol[xx, i]
                ru <- unique(rx)
                if (length(ru) > 1) {
                    ui <- NA
                } else {
                    ui <- ru
                }
                #ui <- ifelse(length(ru) > 1, NA, ru)
                ul <- rep(ui, ll[x])
                return(ul)
            })
            rnCol[, i] <- unlist(ri)

        }
        rownames(rnCol) <- names(idR)
    } else { rnCol <- NULL }


    ## Create the link data
    lvr <- rep(level, ll)
    tipN <- setdiff(ed[, 2], ed[, 1])
    lkr <- DataFrame(nodeLab = transNode(tree = tree,
                                         node = lvr,
                                         use.alias = FALSE,
                                         message = FALSE),
                     nodeLab_alias = transNode(tree = tree,
                                               node = lvr,
                                               use.alias = TRUE,
                                               message = FALSE),
                     nodeNum = lvr,
                     isLeaf = lvr %in% tipN)

    rnD <- LinkDataFrame(LinkData = lkr, rnCol)

    out <- list(dataTab = outR, LinkDF = rnD)
    return(out)
}






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
#' @param rowTree An optional argument. The tree structure on the rows of the
#'   matrix. If required, a \code{phylo} object should be provided.
#' @param colTree An optional argument. The tree structure on the columns of the
#'   matrix. If required, a \code{phylo} object should be provided.
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
#'
#' @importFrom SummarizedExperiment colData rowData 'colData<-' 'rowData<-'
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
#' aggRow <- aggValue(x = tse, colLevel = c("GroupA", "GroupB"),
#' FUN = sum)
#'
#' assays(aggRow)[[1]]
setGeneric("aggValue", function(x, rowTree = NULL, colTree = NULL,
                                rowLevel = NULL, colLevel = NULL,
                                FUN = sum, message = FALSE) {
    standardGeneric("aggValue")
})


#' @describeIn aggValue Aggregation with input x as a matrix
#' @importFrom utils flush.console
setMethod("aggValue", signature(x = "matrix"),
          .aggMatrix)

#' @describeIn aggValue Aggregation with input x as a TreeSummarizedExperiment
#' @importFrom methods is
#' @importFrom utils flush.console
setMethod("aggValue", signature(x = "TreeSummarizedExperiment"),
          .aggTSE)

