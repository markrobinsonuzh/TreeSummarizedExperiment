
#' Perform data aggregations based on the available tree structures
#'
#' \code{aggValue} aggregates values on the leaf nodes of a tree to a specific
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
#' @param assay A integer scalar or string indicating which assay of \code{x} to
#'   use in the aggregation.
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
#'
#' tse <- TreeSummarizedExperiment(assays = list(toyTable),
#'                                 colData = colInf,
#'                                 colTree = treeC,
#'                                 colNodeLab = treeC$tip.label,
#'                                 metadata = list(test = 1:4))
#'
#' aggCol <- aggValue(x = tse, colLevel = c("GroupA", "GroupB"),
#' FUN = sum)
#'
#' assays(aggCol)[[1]]
#'
aggValue <- function(x, rowLevel = NULL, rowBlock = NULL, 
                     colLevel = NULL, colBlock = NULL, 
                     FUN = sum, assay = NULL,
                     message = FALSE) {
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
    if (!is.null(assay)) {
        mat <- mat[assay]
    }
    
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
        if (message) {
            message("The row aggregation is using ", deparse(substitute(FUN)))
        }
        
        rD <- rowData(x)
        rL <- rowLinks(x)
        outR <- .aggFun(tree = rTree,
                        assayTab = mat,
                        dimData = rD,
                        linkData = rL,
                        level = unique(rowLevel),
                        block = rowBlock,
                        FUN = FUN,
                        message = message)
        nrD <- outR$newDD
        nrD <- nrD[, setdiff(colnames(nrD), "nodeLab")]
        nrL <- outR$newDD$nodeLab
        mat <- outR$dataTab
    }else{
        nrD <- rowData(x)
        if (!is.null(rTree)) {
            nrL <- rowLinks(x)$nodeLab
        } else {
            nrL <- NULL
        }
    }
    
    ## -------------------- aggregation on column dimension ----------------
    if (onCol) {
        if (message) {
            message("The column aggregation is using ",
                    deparse(substitute(FUN)))
        }
        cD <- colData(x)
        cL <- colLinks(x)
        outC <- .aggFun(tree = cTree,
                        assayTab = lapply(mat, t),
                        dimData = cD,
                        linkData = cL,
                        level = unique(colLevel),
                        block = colBlock,
                        FUN = FUN, message = message)
        ncD <- outC$newDD
        ncD <- ncD[, setdiff(colnames(ncD), "nodeLab")]
        ncL <- outC$newDD$nodeLab
        mat <- lapply(outC$dataTab, t)
    } else {
        ncD <- colData(x)
        
        if (!is.null(cTree)) {
            ncL <- colLinks(x)$nodeLab
        } else {ncL <- NULL}
    }
    
    # metadata
    metaD <- metadata(x)
    # create the new TreeSummarizedExperiment object
    out <- TreeSummarizedExperiment(assays = mat,
                                    rowData = nrD,
                                    colData = ncD,
                                    rowTree = rTree,
                                    rowNodeLab = nrL,
                                    colTree = cTree,
                                    colNodeLab = ncL,
                                    metadata = metaD)
    
    return(out)
}

########################### version 2 ##############################

.aggFun <- function(tree, assayTab, dimData, linkData,
                    level = NULL, block = NULL, FUN, 
                    message = FALSE) {
    
    # nodeNum & block
    numR <- linkData$nodeNum
    if (!is.null(block)) {
        bk <- dimData[[block]]  
    } else {
        bk <- rep(1, length(numR))
    }
    
    nodeBk <- data.frame(numR, bk, stringsAsFactors = FALSE)
    
    
    # All node numbers on the tree
    ed <- tree$edge
    numAR <- unique(as.vector(ed))
    
    # The aggregation level
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
    miR <- setdiff(unique(unlist(desR)), numR)
    if (length(miR)) {
        warning(length(miR), 
                " leaves couldn't be found from the assay table.\n")
        miR <- transNode(tree = tree, node = miR,
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
    
    assayL <- annL <- vector("list", length(desR))
    ll <- rep(NA, length(desR))
    for (i in seq_along(desR)) {
        if (message) {
            message(i, " out of ", length(desR),
                    " finished", "\r", appendLF = FALSE)
            flush.console()
        }
        
        # assay table
        ri <- which(numR %in% desR[[i]])
        sri <- split(ri, nodeBk[ri, "bk"])
        
        ax <- lapply(sri, FUN = function(x) {
            xx <- mtab[x, , drop = FALSE]
            fx <- apply(xx, 2, FUN = FUN)
            fx <- rbind(fx)
            rownames(fx) <- NULL
            fx
        })
        
        assayL[[i]] <- do.call(rbind, ax)
        
        #dimension data
        li <- lapply(ax, nrow)
        ll[[i]] <- sum(unlist(li))
        
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
            annL[[i]] <- do.call(rbind, ru)
        }}
    
    # separate the results into the original number of assays tables as input
    if (message) {
        message("unwrap data ... ")}
    
    assayL <- do.call(rbind, assayL)
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
    } else {
        newDD <- dimData[rep(1, sum(ll)), , drop = FALSE]
    }
    
    newDD$nodeLab <- rep(names(desR), ll)
    newDD <- DataFrame(newDD)
    rownames(newDD) <- rep(names(desR), ll)
    
    out <- list(dataTab = outR, newDD = newDD)
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
    df[rep(seq_len(nrow(df)), each = n), ]
}
