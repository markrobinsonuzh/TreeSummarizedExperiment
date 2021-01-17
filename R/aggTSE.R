#' Perform data aggregations based on the available tree structures
#'
#' \code{aggTSE} aggregates values on the leaf nodes of a tree to a specific
#' arbitrary level of the tree. The level is specified via the nodes of the
#' tree. Users could decide on which dimension (row or column) and how should
#' the aggregation be performed.
#'
#' @param x A \code{TreeSummarizedExperiment} object.
#' @param rowLevel A numeric (node numbers) or character (node labels) vector.
#'   It provides the level on the tree that data is aggregated to. The
#'   aggregation is on the row dimension. The default is \code{rowLevel = NULL},
#'   and no aggregation is performed.
#' @param rowBlock A column name in the \code{rowData} to separate the
#'   aggregation.
#' @param colLevel A numeric (node numbers) or character (node labels) vector.
#'   It provides the level on the tree that data is aggregated to. The
#'   aggregation is on the column dimension. The default is \code{colLevel =
#'   NULL}, and no aggregation is performed.
#' @param colBlock A column name in the \code{colData} to separate the
#'   aggregation.
#' @param rowFun A function to be applied on the row aggregation. It's similar
#'   to the \code{FUN} in \code{\link[base]{apply}}.
#' @param colFun A function to be applied on the col aggregation. It's similar
#'   to the \code{FUN} in \code{\link[base]{apply}}.
#' @param whichAssay A integer scalar or string indicating which assay of
#'   \code{x} to use in the aggregation. If NULL, all assay tables are used in
#'   aggregation.
#' @param whichRowTree A integer scalar or string indicating which row tree is
#'   used in the aggregation. The first row tree is used as default.
#' @param whichColTree A integer scalar or string indicating which row tree is
#'   used in the aggregation. The first row tree is used as default.
#' @param whichRowTree A integer scalar or string indicating which row tree is
#'   used in the aggregation. The first row tree is used as default.
#' @param rowDataCols The rowData columns to include.
#' @param colDataCols The colData columns to include.
#' @param rowFirst TRUE or FALSE. If the aggregation is in both dims., it is
#'   performed firstly on the row dim for \code{rowFirst = TRUE} or on the
#'   column dim for \code{rowFirst = FALSE}.
#' @param message A logical value. The default is TRUE. If TRUE, it will print
#'   out the running process.
#' @param BPPARAM Default is \code{NULL} and the computation isn't run in
#'   parallel. To run computation parallelly, an optional
#'   \code{\link[BiocParallel:BiocParallelParam-class]{BiocParallelParam}}
#'   instance determining the parallel back-end to be used during evaluation, or
#'   a list of BiocParallelParam instances, to be applied in sequence for nested
#'   calls to \strong{BiocParallel} functions.
#' @import SingleCellExperiment
#' @importFrom utils flush.console
#' @importFrom BiocParallel bplapply
#' @include allClass.R
#' @export
#' @author Ruizhu HUANG
#' @return A \code{\link{TreeSummarizedExperiment}} object 
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
#'
#' tse <- TreeSummarizedExperiment(assays = list(toyTable),
#'                                 colData = colInf,
#'                                 colTree = treeC,
#'                                 colNodeLab = treeC$tip.label,
#'                                 metadata = list(test = 1:4))
#'
#' aggCol <- aggTSE(x = tse, colLevel = c("GroupA", "GroupB"),
#'                  colFun = sum)
#'
#' assays(aggCol)[[1]]
#'
aggTSE <- function(x, 
                   rowLevel = NULL, rowBlock = NULL, 
                   colLevel = NULL, colBlock = NULL, 
                   rowFun = sum, colFun = sum,
                   whichRowTree = 1, whichColTree = 1,
                   whichAssay = NULL,
                   message = FALSE,
                   rowDataCols,
                   colDataCols,
                   rowFirst = TRUE, BPPARAM = NULL) {
    ## --------------------- check input before run ---------------------------
    ## x is a TreeSummarizedExperiment
    if (!is(x, "TreeSummarizedExperiment")) {
        stop("x should be a TreeSummarizedExperiment", call. = FALSE)
    }
    x <- updateObject(x)
    
    ## The provided rowLevel or colLevel should be a numeric or character vector
    isR <- (is.character(rowLevel) | is.numeric(rowLevel) |
                is.null(rowLevel))
    isC <- (is.character(colLevel) | is.numeric(colLevel) |
                is.null(colLevel))
    if (!isR) {
        stop("rowLevel should be a numeric or character vector. \n",
             call. = FALSE)
    }
    if (!isC) {
        stop("colLevel should be a numeric or character vector. \n",
             call. = FALSE)
    }
    
    ## ---------------------- get all data ready ------------------------------
    if (message) {message("Preparing data... ")}
    
    ## The trees
    rTree <- rowTree(x, whichTree = whichRowTree)
    cTree <- colTree(x, whichTree = whichColTree)
    
    ## subset x by trees
    rtn <- rowTreeNames(x)[whichRowTree]
    ctn <- colTreeNames(x)[whichColTree]
    
    ## Indicate whether aggregation should be on rows or columns
    if (is.null(rowLevel)) {
        onRow <- FALSE
    } else {
        if (is.null(rTree)) {stop("The tree on rows doesn't exist.",
                                  call. = FALSE)}
        # subset x by rows that are mapped to the specified tree
        ii <- rowLinks(x)$whichTree == rtn
        x <- x[ii, ]
        si <- sum(ii) < nrow(x)
        if (message & si) {
            message("Subset x by tree: ", sum(ii), " rows are left.")
        }
        onRow <- TRUE
        if (message) {
            message("Perform aggregation on the row dimension... ")
        }
    }
    
    if (is.null(colLevel)) {
        onCol <- FALSE
    } else {
        if (is.null(cTree)) {stop("The tree on columns doesn't exist.",
                                  call. = FALSE)}
        
        # subset x by rows that are mapped to the specified tree
        jj <- colLinks(x)$whichTree == ctn
        x <- x[, jj]
        sj <- sum(jj) < ncol(x)
        if (message & sj) {
            message("Subset x by tree: ", sum(jj), " cols are left.")
        }
        onCol <- TRUE
        if (message) {
            message("Perform aggregation on the column dimension... ")
        }
    }
    
    if (missing(rowDataCols)) {
        rowDataCols <- colnames(rowData(x))
    }
    if (missing(colDataCols)) {
        colDataCols <- colnames(colData(x))
    }
    
    ## The assay tables
    mat <- assays(x)
    if (!is.null(whichAssay)) {mat <- mat[whichAssay]}
    
    ## -------------------- aggregation on row dimension -------------------
    if (rowFirst) {
        x <- .agg_row(x = x, tree = rTree, assayTab = mat,
                      level = unique(rowLevel),
                      rowDataCols = rowDataCols,
                      block = rowBlock, FUN = rowFun, 
                      message = message, onRow = onRow,
                      BPPARAM = BPPARAM)
        ## Update the assay tables for the aggregation on the other dim.
        mat <- assays(x)
    }
    ## -------------------- aggregation on column dimension ----------------
    x <- .agg_col(x = x, tree = cTree, assayTab = mat,
                  level = unique(colLevel),
                  colDataCols = colDataCols,
                  block = colBlock, FUN = colFun, 
                  message = message, onCol = onCol,
                  BPPARAM = BPPARAM)
    
    if (!rowFirst) {
        ## Update the assay tables
        mat <- assays(x)
        x <- .agg_row(x = x, tree = rTree, assayTab = mat,
                      level = unique(rowLevel),
                      rowDataCols =rowDataCols,
                      block = rowBlock, FUN = rowFun, 
                      message = message, onRow = onRow,
                      BPPARAM = BPPARAM)
    }
    
    updateObject(x)
    return(x)
}


.agg_row <- function(x, tree, assayTab, level,
                     rowDataCols,
                     block, FUN, message, onRow, BPPARAM) {
    
    if (onRow) {
        if (message) {
            message("The row aggregation is using ", deparse(substitute(FUN)))
        }
        
        rD <- rowData(x)
        rL <- rowLinks(x)
        outR <- .aggFun(tree = tree,
                        assayTab = assayTab,
                        dimData = rD,
                        linkData = rL,
                        cols = rowDataCols,
                        level = level,
                        block = block,
                        FUN = FUN,
                        message = message,
                        BPPARAM = BPPARAM)
        nrD <- outR$newDD
        nrD <- nrD[, setdiff(colnames(nrD), "nodeLab"), drop = FALSE]
        nLab <- outR$newDD$nodeLab
        mat <- outR$dataTab
    } else {
        mat <- assayTab
        nrD <- rowData(x)[, rowDataCols, drop = FALSE]
        if (!is.null(tree)) {
            nLab <- rowLinks(x)$nodeLab
        } else {
            nLab <- NULL
        }
    }
    
    # new TSE: 
    xx <- TreeSummarizedExperiment(assays = mat,
                                   rowData = nrD,
                                   rowTree = tree,
                                   rowNodeLab = nLab,
                                   referenceSeq = NULL,
                                   colData = colData(x),
                                   metadata = metadata(x))
    xx <- BiocGenerics:::replaceSlots(object = xx, 
                                      colTree = colTree(x, whichTree = NULL),
                                      colLinks = colLinks(x),
                                      check = FALSE)
    return(xx)
}



.agg_col <- function(x, tree, assayTab, level,
                     colDataCols,
                     block, FUN, message, onCol,
                     BPPARAM) {
    
    if (onCol) {
        if (message) {
            message("The column aggregation is using ",
                    deparse(substitute(FUN)))
        }
        cD <- colData(x)
        cL <- colLinks(x)
        outC <- .aggFun(tree = tree,
                        assayTab = lapply(assayTab, t),
                        dimData = cD,
                        linkData = cL,
                        cols = colDataCols,
                        level = level,
                        block = block,
                        FUN = FUN, 
                        message = message,
                        BPPARAM = BPPARAM)
        ncD <- outC$newDD
        ncD <- ncD[, setdiff(colnames(ncD), "nodeLab"), drop = FALSE]
        nLab <- outC$newDD$nodeLab
        mat <- lapply(outC$dataTab, t)
    } else {
        mat <- assayTab
        ncD <- colData(x)[, colDataCols, drop = FALSE]
        if (!is.null(tree)) {
            nLab <- colLinks(x)$nodeLab
        } else {nLab <- NULL}
        
    }
    
    xx <- TreeSummarizedExperiment(assays = mat,
                                   colData = ncD,
                                   colTree = tree,
                                   colNodeLab = nLab,
                                   referenceSeq = referenceSeq(x),
                                   rowData = rowData(x),
                                   metadata = metadata(x))
    xx <- BiocGenerics:::replaceSlots(object = xx, 
                                      rowTree = rowTree(x, whichTree = NULL),
                                      rowLinks = rowLinks(x),
                                      check = FALSE)
    return(xx)
}

.aggFun <- function(tree, assayTab, dimData, linkData,
                    cols, 
                    level = NULL, block = NULL, FUN, 
                    message = FALSE, BPPARAM) {
    
    # nodeNum & block
    numR <- linkData$nodeNum
    if (!is.null(block)) {
        bk <- dimData[[block]]  
    } else {
        bk <- rep(1, length(numR))
    }
    
    nodeBk <- data.frame(numR, bk, stringsAsFactors = FALSE)
    
    # dimData
    dimData <- dimData[, cols, drop = FALSE]
    
    # All node numbers on the tree
    ed <- tree$edge
    numAR <- unique(as.vector(ed))
    
    # The aggregation level
    if (is.character(level)) {
        level <- convertNode(tree = tree, node = level,
                             use.alias = FALSE, message = FALSE)
    }
    
    # The descendants of the aggregation level
    desR <- findDescendant(tree = tree, node = level,
                           only.leaf = TRUE, self.include = TRUE)
    names(desR) <- convertNode(tree = tree, node = level,
                               use.alias = TRUE, message = FALSE)
    
    
    # Find the rows of the descendants
    miR <- setdiff(unique(unlist(desR)), numR)
    if (length(miR)) {
        warning(length(miR), 
                " leaves couldn't be found from the assay table.\n")
        miR <- convertNode(tree = tree, node = miR,
                           use.alias = FALSE, message = FALSE)
        warning("Missing leaves: ", paste(head(miR), collapse = " "), " ...")
    }
    
    ## Perform the aggregation
    # on the assay tables
    
    if (message) {
        message("Working on the assays table... ")}
    
    # aggregate all assay tables together
    nc <- ncol(assayTab[[1]])
    ncn <- colnames(assayTab[[1]])
    mtab <- do.call(cbind, assayTab)
    
    if (is.null(BPPARAM)) {
        # The computation is not run in parallel
        # lapply is in some case much faster than `bplapply` 
        # (https://gist.github.com/fionarhuang/c59c0d1eabd36d8b6fbfda61d06cb951)
        
        res <- lapply(seq_along(desR), FUN = function(x) {
            if (message) {
                message(x, " out of ", length(desR),
                        " finished", "\r", appendLF = FALSE)
                flush.console()
            }
            .agg_assay(desd = desR[[x]], nodeRow = numR, 
                       nodeBlock = nodeBk, assayTab = mtab,
                       dimData = dimData, FUN = FUN)})
    } else {
        # The computation is run in parallel
        res <- bplapply(desR, FUN = function(x) {
            .agg_assay(desd = x, nodeRow = numR, 
                       nodeBlock = nodeBk, assayTab = mtab,
                       dimData = dimData, FUN = FUN)},
            BPPARAM = BPPARAM)
    }
    
    
    # separate the results into the original number of assays tables as input
    if (message) {
        message("unwrap data ... ")}
    
    assayL <- do.call(rbind, lapply(res, FUN = function(x) {x$new_assays}))
    ll <- lapply(res, FUN = function(x) {x$new_n})
    ll <- unlist(ll)
    annL <- lapply(res, FUN = function(x) {x$new_dimData})
    
    outR <- vector("list", length(assayTab))
    names(outR) <- names(assayTab)
    for (i in seq_along(outR)) {
        ii <- seq(from = (i-1)*nc + 1, to = i*nc, by = 1)
        rii <- assayL[, ii, drop = FALSE]
        colnames(rii) <- ncn
        outR[[i]] <- rii
    }
    
    # use the same row names as outR
    if (ncol(dimData)) {
        newDD <- do.call(rbind, annL)
        colnames(newDD) <- cols
    } else {
        newDD <- dimData[rep(1, sum(ll)), , drop = FALSE]
    }
    
    newDD$nodeLab <- rep(names(desR), ll)
    newDD <- DataFrame(newDD)
    rownames(newDD) <- rep(names(desR), ll)
    
    out <- list(dataTab = outR, newDD = newDD)
    return(out)
}

.agg_assay <- function(desd, nodeRow, nodeBlock, assayTab,
                       dimData, FUN) {
    # assay table
    ri <- which(nodeRow %in% desd)
    sri <- split(ri, nodeBlock[ri, "bk"])
    
    ax <- lapply(sri, FUN = function(x) {
        xx <- assayTab[x, , drop = FALSE]
        fx <- apply(xx, 2, FUN = FUN)
        fx <- rbind(fx)
        rownames(fx) <- NULL
        fx
    })
    rax <- do.call(rbind, ax)
    
    
    li <- lapply(ax, nrow)
    sli <- sum(unlist(li))
    
    if (ncol(dimData)) {
        dl <- lapply(sri, FUN = function(x) {
            if (!length(x)) {
                xx <- dimData[1, , drop = FALSE]
                xx[] <- NA
                xx
            } else {
                xx <- as.list(dimData[x, , drop = FALSE])
                ux <- lapply(xx, .unique_or_na)
                do.call(cbind.data.frame, ux)
            }
        })
        ru <- mapply(.rep_df, df = dl, n=li, SIMPLIFY = FALSE)
        dimData <- do.call(rbind, ru)
    }
    
    out <- list(new_assays = rax,
                new_n = sli,
                new_dimData = dimData)
    return(out)
}


.unique_or_na <- function(x) {
    xx <- unique(x)
    if (length(xx) > 1) {
        NA
    } else {
        xx
    }
}

.rep_df <- function(df, n){
    dff <- df[rep(seq_len(nrow(df)), each = n), ]
    data.frame(dff)
}
