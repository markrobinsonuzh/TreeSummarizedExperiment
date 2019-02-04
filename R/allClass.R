
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
#' adding one slot \code{LinkData}
#' @importClassesFrom S4Vectors DataFrame
#' @exportClass LinkDataFrame
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
#' @importFrom S4Vectors DataFrame
#' @export
#' @return A LinkDataFrame object
#'
LinkDataFrame <- function(LinkData = NULL, ...) {
    df <- DataFrame(...)
    new("LinkDataFrame", df, LinkData = LinkData)
}

#-------------------------------------------------------------------------------
# treeSummarizedExperiment
#-------------------------------------------------------------------------------
#' An S4 class treeSummarizedExperiment
#'
#' The class \strong{treeSummarizedExperiment} is an extension class of standard
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} class. It has six
#' slots. Four of them are traditional slots from
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} class:
#' \code{assays}, \code{rowData} \code{colData} and \code{metadata}. The other
#' two slots are \code{linkData} and \code{treeData}. The class
#' \strong{treeSummarizedExperiment} is designed to store rectangular data for
#' entities (e.g., microbes or cell types) (\code{assays}), information about
#' the hiearchical structure of entities (\code{treeData}), and information
#' about the mapping between the rows of the rectangular data and the nodes of
#' the tree (\code{linkData}).
#'
#' @slot linkData A \code{\link[S4Vectors]{DataFrame-class}} object. It gives
#'   map information between the rows of rectangular data and the nodes of tree.
#'   \itemize{
#'   \item \strong{nodeLab} The node labels on the tree.
#'   \item \strong{nodeLab_alias} An alias of column \code{nodeLab}. It is
#'   created only when there are missing value or duplicated value in column
#'   \code{nodeLab}. A prefix "Node_" and "Leaf_" is added to the node number
#'   (column \code{nodeNum}) for the internal nodes and the leaf nodes,
#'   respectively.
#'   \item \strong{nodeNum} The node numbers on the tree.
#'   \item \strong{isLeaf} This indicates whether a node is a leaf node.
#'   \item \strong{rowID} The row number in \code{assays}.
#'   }
#' @slot treeData A phylo object. It gives information about the hiearchical
#'   structure of the entities.
#' @slot ... See \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
#'   for more details about the slots inherited from \code{SummarizedExperiment}
#'   class.
#'
#' @section Constructor:
#' See \code{\link{treeSummarizedExperiment-constructor}} for constructor
#' functions.
#'
#' @section Accessor:
#' See \code{\link{treeSummarizedExperiment-accessor}} for accessor functions.
#'
#' @importFrom methods setClass
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame
#' @name treeSummarizedExperiment-class
#' @exportClass treeSummarizedExperiment
#' @seealso \code{\link{treeSummarizedExperiment}}
#'   \code{\link{treeSummarizedExperiment-accessor}}
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment-class}}
setClass("treeSummarizedExperiment",
         representation(linkData = "DataFrame",
                        treeData = "phylo"),
         contains = "SummarizedExperiment",
         validity = checkTSE)




# ==========================================================================
### Constructor
# ==========================================================================

#' Construct a treeSummarizedExperiment object
#'
#' \code{treeSummarizedExperiment} constructs a treeSummarizedExperiment object.
#'
#' @param tree A phylo object
#' @param linkData A data frame or \code{\link[S4Vectors]{DataFrame-class}}. It
#'   has the same number of rows as the matrix-like elements of \code{assays}.
#'   The row order of the \code{linkData} matches with that of  the matrix-like
#'   element of \code{assays}. It has the following columns.
#'   \itemize{
#'   \item \strong{nodeLab} The labels of nodes on the tree.
#'   \item \strong{nodeLab_alias} An alias of column \code{nodeLab}. It is
#'   created only when there are missing value or duplicated value in column
#'   \code{nodeLab}. A prefix "Node_" and "Leaf_" is added to the node number
#'   (column \code{nodeNum}) for the internal nodes and the leaf nodes,
#'   respectively.
#'   \item \strong{nodeNum} The numbers of nodes
#'   \item \strong{isLeaf} It indicates whether the node is a leaf node
#'   \item \strong{rowID} The corresponding row number of the matrix-like
#'   elements in \code{assays} }
#' @inheritParams SummarizedExperiment::SummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment rowData
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new
#' @export
#' @include classValid.R
#' @return a treeSummarizedExperiment object
#' @name treeSummarizedExperiment-constructor
#' @author Ruizhu HUANG
#' @seealso \code{\link{treeSummarizedExperiment-class}}
#'   \code{\link{treeSummarizedExperiment-accessor}}
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
#' # build a treeSummarizedExperiment object
#' tse <- treeSummarizedExperiment(tree = tinyTree,
#'                                 assays = list(count),
#'                                 colData = sampC)
#'
treeSummarizedExperiment <- function(tree = NULL, linkData = NULL,
                                     ...){

    # -------------------------------------------------------------------------
    # create a SummarizedExperiment object
    se <- SummarizedExperiment(...)

    # -------------------------------------------------------------------------
    # tree is a phylo object
    if(!inherits(tree, "phylo")) {
        stop(tree, ": A phylo object is required", "\n")
    }

    # -------------------------------------------------------------------------
    # The labels of tree nodes should be unique
    treeLab <- c(tree$tip.label, tree$node.label)
    tipLab <- tree$tip.label
    isDp <- duplicated(tipLab)
    anyDp <- any(isDp)
    if (anyDp) {
        stop("\n Can not distinguish nodes with the same label: ",
             head(tipLab[isDp])," \n")
    }

    # -------------------------------------------------------------------------
    ### create link data
    if (is.null(linkData)) {

        # (1) if the nodeLab column exist, they should match with the labels of
        # tree nodes.
        nodeLab <- rowData(se)$nodeLab
        if (!is.null(nodeLab)) {
            # keep only rows that could be assigned to the nodes of tree
            isIn <- nodeLab %in% treeLab
            isOut <- !isIn
            if (sum(isOut) > 0) {
                message(sum(isOut), "rows are removed. \n",
                        "They can't be matched to any node of the tree. \n")}
            se <- se[isIn, ]
            newLab <- rowData(se)$nodeLab
            }

        # (2) if the nodeLab column doesn't exist, rownames should match with
        # the labels of tree leaves.
        if (is.null(nodeLab)) {
            rowNam <- rownames(se)

            # if neither nodeLab nor rownames are provided, return error.
            if (is.null(rowNam)) {
                stop("Either a nodeLab column or row names should be
                     provided for row data \n.")}
            # keep only rows that could be assigned to the nodes of tree
            isIn <- rowNam %in% treeLab
            isOut <- !isIn
            if (sum(isOut) > 0) {
                message(sum(isOut), " rows are removed. They cannot be maped to
                    any node of the tree. \n")}
            se <- se[isIn, ]
            newLab <- rownames(se)
            }

        # create linkD
        linkD <- DataFrame(nodeLab = newLab,
                           nodeNum = transNode(tree = tree,
                                               input = newLab,
                                               message = FALSE),
                           isLeaf = newLab %in% tree$tip.label,
                           rowID = seq_len(length(newLab)))

        # create column nodeLab_alias, if there are duplicated value in the
        # nodeLab column
        if (any(duplicated(linkD$nodeLab))) {
            linkD$nodeLab_alias <- transNode(tree = tree,
                                             input = linkD$nodeNum,
                                             use.alias = TRUE)
        }

        } else {
            # if linkData is provided, then use it as linkData
            linkD <- linkData

            # create column nodeLab_alias, if there are duplicated value in the
            # nodeLab column
            if (any(duplicated(linkD$nodeLab))) {
                linkD$nodeLab_alias <- transNode(tree = tree,
                                                 input = linkD$nodeNum,
                                                 use.alias = TRUE)
            }
        }

    # -------------------------------------------------------------------------
    # create treeSummarizedExperiment
    rowData(se) <- rowData(se)[, colnames(rowData(se)) != "nodeLab",
                               drop = FALSE]
    tse <- new("treeSummarizedExperiment", se,
               linkData = linkD, treeData = tree)

    return(tse)
}

