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
        if (!inherits(object@treeData, "phylo")) {
            msg <- cat("\n tree is not a phylo object")
        }

        # The leaf nodes should have unique label.
        tree <- object@treeData
        tipLab <- tree$tip.label
        isDp <- duplicated(tipLab)
        anyDp <- any(isDp)
        if (anyDp) {
            msg <- cat("\n Different leaf nodes using the same label: ",
                       head(tipLab[isDp])," \n")
            errors <- c(errors, msg)
        }
    }


    # -------------------------------------------------------------------------
    # if the linkData exists, a column nodeLab_alias should be generated if
    # there are duplicated value in the nodeLab column
    lk <- object@linkData
    if (!is.null(lk)) {
        if (anyDuplicated(lk$nodeLab)) {
            if (is.null(lk$nodeLab_alias)) {
                msg <- cat("\n Duplicated values in the column nodeLab. \n")
                errors <- c(errors, msg)
            }
        }
    }

    # -------------------------------------------------------------------------
    # Note : duplicated value in nodeLab column is allowed because we might
    # have multiple rows corresponding to a same leaf.
    if (length(errors) == 0) {TRUE} else {errors}
}
