#' valid TreeSummarizedExperiment class
#'
#' \code{validTSE} is to valid TreeSummarizedExperiment object.
#'
#'
#' @param object A TreeSummarizedExperiment object
#'
#' @import SingleCellExperiment S4Vectors
#' @return TRUE or a character string.
#' @keywords internal
#'
.checkTSE <- function(object){

    errors <- character()
    # -------------------------------------------------------------------------
    # If provided, the row tree should be a phylo object; otherwise, NULL
    rTree <- rowTree(object)

    if (!is.null(rTree)) {
        if (!is(rTree, "phylo")) {
        msg <- "The row tree is not a phylo object. \n"
        errors <- c(errors, msg)
        }
     }

    # If provided, the column tree should be a phylo object; otherwise, NULL
    cTree <- colTree(object)
    if (!is.null(cTree)) {
        if (!is(cTree, "phylo")) {
            msg <- "The column tree is not a phylo object. \n"
            errors <- c(errors, msg)
        }
    }

    # -------------------------------------------------------------------------
    # The labels of the tree leaves should be unique
    if (!is.null(rTree)) {
        # The leaf nodes should have unique label.
        tipLab <- rTree$tip.label
        anyDp <- any(duplicated(tipLab))

        if (anyDp) {
            msg <- "rowTree: Duplicated labels are not allowed for leaves. \n"
            errors <- c(errors, msg)
        }
    }

    if (!is.null(cTree)) {
        # The leaf nodes should have unique label.
        tipLab <- cTree$tip.label
        anyDp <- any(duplicated(tipLab))

        if (anyDp) {
            msg <- "colTree: Duplicated labels are not allowed for leaves. \n"
            errors <- c(errors, msg)
        }
    }

    # -------------------------------------------------------------------------
    # check the dimension is correct for rowLinks and colLinks
    if (!is.null(rTree)) {
        rowEq <- nrow(object@rowLinks) == nrow(object)

        if (!rowEq) {
            msg <- sprintf("rowLinks: %d rows are expected",
                           nrow(object))
            errors <- c(errors, msg)
        }
        if (!all(rownames(object@rowLinks) == rownames(object))) {
            msg <- "rowLinks: rownames do not match rownames of experiment"
            errors <- c(errors, msg)
        }
    }

    if (!is.null(cTree)) {
        colEq <- nrow(object@colLinks) == ncol(object)

        if (!colEq) {
            msg <- sprintf("colLinks: %d cols are expected",
                           ncol(object))
            errors <- c(errors, msg)
        }
        if (!all(rownames(object@colLinks) == colnames(object))) {
            msg <- "colLinks: rownames do not match colnames of experiment"
            errors <- c(errors, msg)
        }
    }

    # -------------------------------------------------------------------------
    # if rowTree doesn't exist, rowLinks should have 0 rows
    if (is.null(rTree)) {
        if (!is.null(object@rowLinks)) {
            msg <- "rowLinks should be NULL when rowTree doesn't exist \n"
            errors <- c(errors, msg)
        }
    }

    if (is.null(cTree)) {
        if (!is.null(object@colLinks)) {
            msg <- "colLinks should be NULL when colTree doesn't exist \n"
            errors <- c(errors, msg)
        }
    }

    # -------------------------------------------------------------------------
    # Note : duplicated value in nodeLab column is allowed because we might
    # have multiple rows corresponding to a same leaf.
    if (length(errors)) {stop("\n", errors)} else {TRUE}
}

# =============================================================================
# valid the LinkDataFrame class
# =============================================================================
#' \code{validLDF} is to valid a \code{LinkDataFrame} object.
#'
#' @param object A \code{LinkDataFrame} object
#' @import SingleCellExperiment S4Vectors
#' @return TRUE or a character string.
#' @keywords internal
#'
.checkLDF <- function(object){

    errors <- character()
    # -------------------------------------------------------------------------
    # it must be a subclass of DataFrame
    if (!is(object, "DataFrame")) {
        msg <- "The object is not a subclass of DataFrame \n"
        errors <- c(errors, msg)
    }

    # -------------------------------------------------------------------------
    # it should at least include nodeLab, nodeLab_alias, nodeNum, isLeaf
    colNam <- colnames(object)
    rqNam <- c("nodeLab", "nodeLab_alias", "nodeNum", "isLeaf")
    isNam <- all(rqNam %in% colNam)
    if (!isNam) {
        msg <- "The object should include at least 4 columns:
                nodeLab, nodeLab_alias, nodeNum, isLeaf \n"
        errors <- c(errors, msg)
    }

    # -------------------------------------------------------------------------
    # nodeLab: a character column (allows NA)
    nodeLab <- object$nodeLab
    isC1 <- is(nodeLab, "character") | all(is.na(nodeLab))
    if (!isC1) {
        msg <- "The nodeLab column should be character \n"
        errors <- c(errors, msg)
    }

    # nodeLab_alias: a character column
    nodeLab_alias <- object$nodeLab_alias
    isC2 <- is(nodeLab_alias, "character")
    if (!isC2) {
        msg <- "The nodeLab_alias column should be character \n"
        errors <- c(errors, msg)
    }

    # nodeNum: a numeric column
    nodeNum <- object$nodeNum
    isC3 <- is(nodeNum, "numeric")
    if (!isC3) {
        msg <- "The nodeNum column should be numeric \n"
        errors <- c(errors, msg)
    }

    # isLeaf: a logical column
    isLeaf <- object$isLeaf
    isC4 <- is(isLeaf, "logical")
    if (!isC4) {
        msg <- "The isLeaf column should be logical \n"
        errors <- c(errors, msg)
    }

    if (length(errors)) {
        stop("\n", errors)
    } else { TRUE }

}
