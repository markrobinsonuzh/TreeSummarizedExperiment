
# ======================
### class
# ======================

#-------------------------------------------------------------------------------
#' phylo: A S3 class for the phylogenetic tree
#-------------------------------------------------------------------------------
#'
#' The \pkg{ape} package does not export its phylo class, probably because it
#' is not really defined formally anywhere. Technically, it is an S3 class
#' extended from the class list. Any exported definitions from the \pkg{ape}
#' package would be preferred to use if available.
#' @importFrom methods setOldClass
#' @keywords internal
phylo <- structure(list(), class = "phylo")
setOldClass("phylo")


#-------------------------------------------------------------------------------
#' LinkDataFrame: A S4 class extended from DataFrame
#-------------------------------------------------------------------------------
#' An S4 class LinkDataFrame
#'
#' The \strong{LinkDataFrame} is extended from the class \strong{DataFrame} by
#' adding one new slot \code{LinkData}
#'
#' @slot LinkData A \link[S4Vectors]{DataFrame-class}. It will be shown in the
#'   left side of the vertical line when print out the generated
#'   \code{LinkDataFrame}.
#' @slot ... Other slots from \link[S4Vectors]{DataFrame-class}.
#'
#' @importClassesFrom S4Vectors DataFrame
#' @exportClass LinkDataFrame
#'
#' @section Constructor:
#' See \code{\link{LinkDataFrame-constructor}} for constructor
#' functions.
#'
#' @section Accessor:
#' See \code{\link{LinkDataFrame-accessor}} for accessor functions.
#'
#' @seealso \code{\link{LinkDataFrame-accessor}}
#'   \code{\link{LinkDataFrame-constructor}} \link[S4Vectors]{DataFrame-class}
#' @name LinkDataFrame-class
setClass("LinkDataFrame",
         representation(LinkData = "DataFrame"),
         contains = "DataFrame")

#-------------------------------------------------------------------------------
#' Construct a LinkDataFrame
#-------------------------------------------------------------------------------
#' Construct a LinkDataFrame object
#'
#' @param LinkData A \link[S4Vectors]{DataFrame-class}.
#' @inheritParams SummarizedExperiment::SummarizedExperiment
#'
#' @importFrom S4Vectors DataFrame
#' @rdname LinkDataFrame-constructor
#' @export
#' @return A LinkDataFrame object
#' @seealso \code{\link{LinkDataFrame-accessor}}
#'   \code{\link{LinkDataFrame-class}} \code{\link[S4Vectors]{DataFrame-class}}
#' @examples
#'
#' left <- DataFrame(left1 = 1:5, left2 = letters[1:5])
#' right <- DataFrame(right1 = sample(letters[1:3], 5, replace = TRUE),
#'                   right2 = sample(c(TRUE, FALSE), 5, replace = TRUE),
#'                   right3 = 11:15)
#'
#' (ld <- LinkDataFrame(LinkData = left, right))
LinkDataFrame <- function(LinkData = NULL, ...) {
    df <- DataFrame(...)
    new("LinkDataFrame", df, LinkData = LinkData)
}

#-------------------------------------------------------------------------------
# TreeSummarizedExperiment
#-------------------------------------------------------------------------------
#' An S4 class TreeSummarizedExperiment
#'
#' The class \strong{TreeSummarizedExperiment} is an extension class of standard
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} class. It has
#' five slots. Four of them are traditional slots from
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} class:
#' \code{assays}, \code{rowData} \code{colData} and \code{metadata}. The new
#' slot is \code{treeData} that is to store the hiearchical information of rows
#' (or columns or both) of \code{assays} tables.
#'
#' @slot treeData A list of phylo objects. It gives information about the hiearchical
#'   structure of rows or columns of \code{assay} tables.
#' @slot ... Other slots from
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#'
#' @section Constructor:
#' See \code{\link{TreeSummarizedExperiment-constructor}} for constructor
#' functions.
#'
#' @section Accessor:
#' See \code{\link{TreeSummarizedExperiment-accessor}} for accessor functions.
#'
#' @details The class \strong{TreeSummarizedExperiment} is designed to store
#'   rectangular data for entities (e.g., microbes or cell types)
#'   (\code{assays}), information about the hiearchical structure
#'   (\code{treeData}), and the mapping information between the rows (or
#'   columns, or both) of the rectangular data and the tree nodes
#'   (\code{linkData}). Users could provide hiearchical information on the rows
#'   or columns (or both) of the \code{assay} tables, and the \code{linkData}
#'   will be automatically integrated as one part of the \code{rowData} or
#'   \code{colData} or both, respectively. When the \code{linkData} is added to
#'   \code{rowData} or \code{colData}, a class \code{LinkDataFrame} is used to
#'   store data instead of \code{DataFrame}. Please see the page
#'   \code{\link{LinkDataFrame}} for more details.
#'
#' @importFrom methods setClass
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame
#' @name TreeSummarizedExperiment-class
#' @exportClass TreeSummarizedExperiment
#' @seealso \code{\link{TreeSummarizedExperiment}}
#'   \code{\link{TreeSummarizedExperiment-accessor}}
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
setClass("TreeSummarizedExperiment",
         representation(treeData = "list"),
         contains = "SummarizedExperiment",
         validity = checkTSE)




# ==========================================================================
### Constructor
# ==========================================================================

#' Construct a TreeSummarizedExperiment object
#'
#' \code{TreeSummarizedExperiment} constructs a TreeSummarizedExperiment object.
#'
#' @param rowTree A phylo object that provides hiearchical information of rows
#'   of assay tables.
#' @param colTree A phylo object that provides hiearchical information of
#'   columns of assay tables.
#' @inheritParams SummarizedExperiment::SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new is
#' @export
#' @include classValid.R
#' @return a TreeSummarizedExperiment object
#' @name TreeSummarizedExperiment-constructor
#' @author Ruizhu HUANG
#' @seealso \code{\link{TreeSummarizedExperiment-class}}
#'   \code{\link{TreeSummarizedExperiment-accessor}}
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#' @examples
#' data("tinyTree")
#'
#' # the count table
#' count <- matrix(rpois(100, 50), nrow = 10)
#' rownames(count) <- c(tinyTree$tip.label)
#' colnames(count) <- paste("C_", 1:10, sep = "_")
#'
#' # The sample information
#' sampC <- data.frame(condition = rep(c("control", "trt"), each = 5),
#'                     gender = sample(x = 1:2, size = 10, replace = TRUE))
#' rownames(sampC) <- colnames(count)
#'
#' # build a TreeSummarizedExperiment object
#' tse <- TreeSummarizedExperiment(rowTree = tinyTree,
#'                                 assays = list(count),
#'                                 colData = sampC)
#'
TreeSummarizedExperiment <- function(rowTree = NULL, colTree = NULL,
                                     ...){

    # -------------------------------------------------------------------------
    ## create a SummarizedExperiment object
    se <- SummarizedExperiment(...)

    # -------------------------------------------------------------------------
    ## Indicate whether hiearchical information is available for rows and
    ## columns of assay tables
    isRow <- !is.null(rowTree)
    isCol <- !is.null(colTree)

    # -------------------------------------------------------------------------
    ## check whether the input tree has the correct form
    if (isRow) {
        # the tree should be a phylo object
        if(!is(rowTree, "phylo")) {
            stop("A phylo object is required for the rowTree", "\n")
        }

        # the tree should have unique leaf labels
        tipLab <- rowTree$tip.label
        anyDp <- any(duplicated(tipLab))
        if (anyDp) {
            stop("rowTree should have unique leaf labels. \n")
        }
    }

    if (isCol) {
        # the tree should be a phylo object
        if(!is(colTree, "phylo")) {
            stop("A phylo object is required for the colTree", "\n")
        }

        # the tree should have unique leaf labels
        tipLab <- colTree$tip.label
        anyDp <- any(duplicated(tipLab))
        if (anyDp) {
            stop("colTree should have unique leaf labels. \n")
        }
    }

    # -------------------------------------------------------------------------
    ## create the link data
    # the rows:
    if (isRow) {
        rowLab <- rowData(se)$nodeLab
        treeLab <- c(rowTree$tip.label, rowTree$node.label)

        # (1) The rowLab should match with the node labels of the rowTree
        if (!is.null(rowLab)) {
            isIn <- rowLab %in% treeLab
            isOut <- !isIn
            if (sum(isOut) > 0) {
                message(sum(isOut), "rows couldn't be matched to the tree and are removed. \n")}
            se <- se[isIn, ]
            rowLab <- rowData(se)$nodeLab
        }

        # (2) If the nodeLab doesn't exist, the rownames should be used instead
        if (is.null(rowLab)) {
            rowLab <- rownames(se)

            if (is.null(rowLab)) {
                stop("Either a nodeLab column or row names should be
                     provided for row data \n.")
            }

            isIn <- rowLab %in% treeLab
            isOut <- !isIn
            if (sum(isOut) > 0) {
                message(sum(isOut), "rows couldn't be matched to the tree and are removed. \n")}
            se <- se[isIn, ]
            rowLab <- rownames(se)
        }

        # create the link data
        rowLink <- DataFrame(nodeLab = rowLab,
                           nodeNum = transNode(tree = rowTree,
                                               input = rowLab,
                                               message = FALSE),
                           isLeaf = rowLab %in% rowTree$tip.label)

        # the column nodeLab_alias is created if a node label is used for
        # different node numbers
        isDn <- duplicated(rowLink$nodeNum)
        isDl <- duplicated(rowLink$nodeLab)
        if (any(isDn != isDl)) {
            rowLink$nodeLab_alias <- transNode(tree = rowTree,
                                             input = rowLink$nodeNum,
                                             use.alias = TRUE)}

        # create the rowData
        rd <- rowData(se)
        rd <- rd[, colnames(rd) != "nodeLab", drop = FALSE]
        rowData(se) <- LinkDataFrame(LinkData = rowLink, rd)

    } else {
        rd <- rowData(se)
        rowData(se) <- rd[, colnames(rd) != "nodeLab", drop = FALSE]

        }


    # the columns:
    if (isCol) {
        colLab <- colData(se)$nodeLab
        treeLab <- c(colTree$tip.label, colTree$node.label)

        # (1) The colLab should match with the node labels of the colTree
        if (!is.null(colLab)) {
            isIn <- colLab %in% treeLab
            isOut <- !isIn
            if (sum(isOut) > 0) {
                message(sum(isOut), "columns couldn't be matched to the tree and are removed. \n")}
            se <- se[, isIn ]
            colLab <- colData(se)$nodeLab
        }

        # (2) If the nodeLab doesn't exist, the rownames should be used instead
        if (is.null(colLab)) {
            colLab <- colnames(se)

            if (is.null(colLab)) {
                stop("Either a nodeLab column or row names should be
                     provided for column data \n.")
            }

            isIn <- colLab %in% treeLab
            isOut <- !isIn
            if (sum(isOut) > 0) {
                message(sum(isOut), "columns couldn't be matched to the tree and are removed. \n")}
            se <- se[, isIn]
            colLab <- colnames(se)
        }

        # create the link data
        colLink <- DataFrame(nodeLab = colLab,
                             nodeNum = transNode(tree = colTree,
                                                 input = colLab,
                                                 message = FALSE),
                             isLeaf = colLab %in% colTree$tip.label)

        # the column nodeLab_alias is created if a node label is used for
        # different node numbers
        isDn <- duplicated(colLink$nodeNum)
        isDl <- duplicated(colLink$nodeLab)
        if (any(isDn != isDl)) {
            colLink$nodeLab_alias <- transNode(tree = colTree,
                                               input = colLink$nodeNum,
                                               use.alias = TRUE)}

        # create the colData
        cd <- colData(se)
        cd <- cd[, colnames(cd) != "nodeLab", drop = FALSE]
        colData(se) <- LinkDataFrame(LinkData = colLink, cd)

    } else {
        cd <- colData(se)
        colData(se) <- cd[, colnames(cd) != "nodeLab", drop = FALSE]

    }


    # -------------------------------------------------------------------------
    # create TreeSummarizedExperiment
    tse <- new("TreeSummarizedExperiment", se,
               treeData = list(rowTree = rowTree,
                               colTree = colTree))

    return(tse)
}

