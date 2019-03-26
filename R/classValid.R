#' valid TreeSummarizedExperiment class
#'
#' \code{validTSE} is to valid TreeSummarizedExperiment object.
#'
#'
#' @param object A TreeSummarizedExperiment object
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData
#' @importFrom S4Vectors metadata
#' @importFrom utils head
#' @return TRUE or a character string.
#' @keywords internal
#'
checkTSE <- function(object){

    errors <- character()
    # -------------------------------------------------------------------------
    # it must have table in assays
    if (length(assays(object)) < 0) {
      msg <- cat("\n No table is available in assays. \n")
      errors <- c(errors, msg)
    }

    # -------------------------------------------------------------------------
    #  Tree should be a phylo object
    if (!is.null(object@treeData)) {
        treeD <- object@treeData
        classD <- lapply(treeD, class)
        classD <- unique(unlist(classD))
        isP <- all(classD %in% c("phylo", "NULL"))
        if (!isP) {
            msg <- cat("\n treeData should be a list of phylo objects")
        }

        # The leaf nodes should have unique label.
        tipLab <- lapply(treeD, FUN = function(x) {x$tip.label})
        isDp <- lapply(tipLab, FUN = function(x){any(duplicated(x))})
        anyDp <- any(unlist(isDp))

        if (anyDp) {
            msg <- cat("\n Duplicated labels are found in the tree. \n")
            errors <- c(errors, msg)
        }
    }

    # -------------------------------------------------------------------------
    # Note : duplicated value in nodeLab column is allowed because we might
    # have multiple rows corresponding to a same leaf.
    if (length(errors) == 0) {TRUE} else {errors}
}
