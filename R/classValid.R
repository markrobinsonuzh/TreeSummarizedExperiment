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
    # it must have table in assays
    if (length(assays(object)) < 0) {
      msg <- "No table is available in assays. \n"
      errors <- c(errors, msg)
    }

    # -------------------------------------------------------------------------
    # If provided, the row tree should be a phylo object; otherwise, NULL
    rTree <- object@rowTree
    isRT <- class(rTree) %in% c("phylo", "NULL")
    if (!isRT) {
        msg <- "The row tree is not a phylo object. \n"
        errors <- c(errors, msg)
    }

    # If provided, the column tree should be a phylo object; otherwise, NULL
    cTree <- object@colTree
    isCT <- class(cTree) %in% c("phylo", "NULL")
    if (!isCT) {
        msg <- "The column tree is not a phylo object. \n"
        errors <- c(errors, msg)
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
    # check the dimension is correct for rowLink and colLink
    if (!is.null(rTree)) {
        rowEq <- nrow(object@rowLink) == nrow(object)

        if (!rowEq) {
            msg <- sprintf("rowLink: %d rows are expected",
                           nrow(object))
            errors <- c(errors, msg)
        }
    }

    if (!is.null(cTree)) {
        colEq <- nrow(object@colLink) == ncol(object)

        if (!colEq) {
            msg <- sprintf("rowLink: %d rows are expected",
                           ncol(object))
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
