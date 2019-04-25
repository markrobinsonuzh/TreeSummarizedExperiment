
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

setClassUnion("listOrNULL", c("list", "NULL"))
#-------------------------------------------------------------------------------
#' LinkDataFrame: A S4 class extended from DataFrame
#-------------------------------------------------------------------------------
#' An S4 class LinkDataFrame
#'
#' The \strong{LinkDataFrame} is extended from the class \strong{DataFrame} to
#' include at least four columns \code{nodeLab}, \code{nodeLab_alias},
#' \code{nodeNum}, and \code{isLeaf}.
#'
#' @importClassesFrom S4Vectors DataFrame
#' @exportClass LinkDataFrame
#'
#' @section Constructor:
#' See \code{\link{LinkDataFrame-constructor}} for constructor
#' functions.
#'
setClass("LinkDataFrame",
         contains = "DataFrame",
         validity = .checkLDF)

#-------------------------------------------------------------------------------
#' Construct a LinkDataFrame
#-------------------------------------------------------------------------------
#' Construct a LinkDataFrame object
#' @param nodeLab A character vector
#' @param nodeLab_alias A character vector
#' @param nodeNum A numeric vector
#' @param isLeaf A logical vector
#' @param ... All arguments accepted by
#'   \code{\link[S4Vectors]{DataFrame-class}}.
#'
#' @importFrom S4Vectors DataFrame
#' @name LinkDataFrame-constructor
#' @export
#' @return A LinkDataFrame object
#' @seealso \code{\link{LinkDataFrame-class}}
#'   \code{\link[S4Vectors]{DataFrame-class}}
#' @examples
#'
#'
#' (ld <- LinkDataFrame(nodeLab = letters[1:5],
#'                      nodeLab_alias = LETTERS[1:5],
#'                      nodeNum = 1:5,
#'                      isLeaf = TRUE,
#'                      right = 1:5))
LinkDataFrame <- function(nodeLab, nodeLab_alias, nodeNum,
                          isLeaf, ...) {
    df <- DataFrame(nodeLab, nodeLab_alias, nodeNum,
                    isLeaf, ...)

    new("LinkDataFrame", df)
}

#-------------------------------------------------------------------------------
# TreeSummarizedExperiment
#-------------------------------------------------------------------------------
#' An S4 class TreeSummarizedExperiment
#'
#' The class \strong{TreeSummarizedExperiment} is an extension class of standard
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} class. It has
#' four more slots that are not in
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} class:
#' \code{rowTree}, \code{rowLinks} \code{colTree} and \code{colLinks}. The
#' hierarchical information of rows (columns) is stored in \code{rowTree}
#' (\code{colTree}) and the link between the rows (columns) of \code{assays}
#' tables and nodes of the tree is given in \code{rowLinks} (\code{colLinks}).
#'
#' @slot rowTree A phylo object or NULL. It gives information about the
#'   hiearchical structure of rows of \code{assays} tables.
#' @slot colTree A phylo object or NULL. It gives information about the
#'   hiearchical structure of columns of \code{assays} tables.
#' @slot rowLinks A LinkDataFrame. It gives information about the link between
#'   the nodes of the \code{rowTree} and the rows of \code{assays} tables.
#' @slot colLinks A LinkDataFrame. It gives information about the link between
#'   the nodes of the \code{colTree} and the columns of \code{assays} tables.
#' @slot ... Other slots from
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}}
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
#'   (\code{rowTree} on rows; \code{colTree} on columns), and the mapping
#'   information between the tree nodes and the rows or the columns of the
#'   rectangular data. Users could provide the hiearchical structure of the
#'   rows, columns or both) of the \code{assays} tables, and the link data will
#'   be automatically generated in \code{rowLinks}, \code{colData} or both,
#'   respectively. It's required that the object in \code{rowLinks} or
#'   \code{colLinks} has the \code{LinkDataFrame} class. Please see the page
#'   \code{\link{LinkDataFrame}} for more details.
#'
#' @importFrom methods setClass
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom S4Vectors DataFrame
#' @name TreeSummarizedExperiment-class
#' @exportClass TreeSummarizedExperiment
#' @seealso \code{\link{TreeSummarizedExperiment}}
#'   \code{\link{TreeSummarizedExperiment-accessor}}
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}}
setClass("TreeSummarizedExperiment",
         representation(rowTree = "listOrNULL",
                        colTree = "listOrNULL",
                        rowLinks = "LinkDataFrame",
                        colLinks = "LinkDataFrame"),
         contains = "SingleCellExperiment",
         validity = .checkTSE)




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
#' @inheritParams SingleCellExperiment::SingleCellExperiment
#'
#' @details The output TreeSummarizedExperiment object has very similar
#'   structure as the
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}}. The
#'   differences are summarized be as below.
#'   \itemize{
#'   \item \strong{rowTree} A slot exists in \code{TreeSummarizedExperiment}
#'   but not in \code{SingleCellExperiment}. It stores the tree structure(s)
#'   that provide(s) hierarchical information of \code{assays} rows or columns
#'   or both.
#'   \item \strong{rowData} If a \code{phylo} object is available in the slot
#'   \code{treeData} to provide the hiearchical information about the rows of
#'   the \code{assays} table, the \code{rowData} would be a
#'   \code{\link{LinkDataFrame-class}} instead of
#'   \code{\link[S4Vectors]{DataFrame-class}}. The data on the right side of the
#'   vertical line provides the link information between the \code{assays} rows
#'   and the tree \code{phylo} object, and could be accessed via
#'   \code{linkData}; The data on the left side is the original \code{rowData}
#'   like \code{SingleCellExperiment} object.
#'   \item \strong{colData} Similar to the explanaition for \strong{rowData} as
#'   above.
#'  }
#'  More details about the \code{LinkDataFrame} in the \code{rowData} or
#'  \code{colData}.
#'  \itemize{
#'  \item nodeLab The labels of nodes on the tree.
#'  \item nodeLab\_alias The alias of node labels on the tree.
#'  \item nodeNum The numbers of nodes on the tree.
#'  \item isLeaf It indicates whether the node is a leaf node or internal node.
#'  }
#'
#' @importFrom SummarizedExperiment assays colData<- rowData<-
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new is as
#' @export
#' @include classValid.R
#' @return a TreeSummarizedExperiment object
#' @name TreeSummarizedExperiment-constructor
#' @author Ruizhu HUANG
#' @seealso \code{\link{TreeSummarizedExperiment-class}}
#'   \code{\link{TreeSummarizedExperiment-accessor}}
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' @examples
#'
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
    sce <- SingleCellExperiment(...)

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
    tse <- new("TreeSummarizedExperiment", sce)
    tse@rowTree <- list(phylo = rowTree)

    # the rows:
    if (isRow) {
        rowD <- rowData(sce)
        tse@rowLinks <- .linkFun(tree = rowTree, sce = sce, onRow = TRUE)
        rowData(tse) <- rowD[, setdiff(colnames(rowD), "nodeLab")]
    } else {
        tse@rowLinks <- LinkDataFrame(nodeLab = character(0),
                                      nodeLab_alias = character(0),
                                      nodeNum = numeric(0),
                                      isLeaf = logical(0))
    }

    tse@colTree <- list(phylo = colTree)
    # the columns:
    if (isCol) {
        colD <- colData(sce)
        tse@colLinks <- .linkFun(tree = colTree, sce = sce, onRow = FALSE)
        colData(tse) <- colD[, setdiff(colnames(colD), "nodeLab")]

    } else {
        tse@colLinks <- LinkDataFrame(nodeLab = character(0),
                                      nodeLab_alias = character(0),
                                      nodeNum = numeric(0),
                                      isLeaf = logical(0))
    }

    return(tse)
}

# An internal function to create the link data and added to the rowData or
# colData
.linkFun <- function(tree, sce, onRow = TRUE){

    if (onRow) {
        annDat <- rowData(sce)
    } else {
        annDat <- colData(sce)
    }

    kw <- ifelse(onRow, "row", "column")

    # ------------------ labels from the tree ---------------------------------
    # the labels and the alias of the labels
    treeLab <- c(tree$tip.label, tree$node.label)
    nodeA <- unique(as.vector(tree$edge))
    treeLab_alias <- transNode(tree = tree, node = nodeA,
                               use.alias = TRUE, message = FALSE)



    # ------------------ labels from the assays table -------------------------
    # The order to be used: nodeLab_alias > nodeLab > the row names
    # nodeLab_alias
    lab <- as.data.frame(annDat)$nodeLab_alias

    # if nodeLab_alias is not available, use nodeLab
    if (is.null(lab)) {
        lab <- as.data.frame(annDat)$nodeLab
    }

    # if nodeLab_alias and nodeLab are not available, use the row names
    if (is.null(lab)) {
        lab <- rownames(annDat)
        if (is.null(lab)) {
            stop("Either a nodeLab column or the row names should be
                 provided for ", kw ," data. \n")}
        }

    # decide whether treeLab or treeLab_alias should be used
    sw <- startsWith(lab, "alias_")
    sw <- all(sw)

    # Match lab with the alias of the node labels on the tree
    if (sw) {
        isIn <- lab %in% treeLab_alias
    } else {
        isIn <- lab %in% treeLab
    }

    # Those with labels that don't match with the node labels of the tree are
    # excluded
    isOut <- !isIn
    if (sum(isOut) > 0) {
        message(sum(isOut), " ", kw,
                "(s) couldn't be matched to the tree and are/is removed. \n")}

    if (onRow) {
        sce <- sce[isIn, ]
    } else {
        sce <- sce[, isIn]
    }


    # create the link data
    labN <- lab[isIn]
    nd <- transNode(tree = tree, node = labN, use.alias = FALSE,
                    message = FALSE)
    fLab <- transNode(tree = tree, node = nd, use.alias = FALSE,
                      message = FALSE)
    faLab <- transNode(tree = tree, node = nd, use.alias = TRUE,
                       message = FALSE)
    leaf <- unique(setdiff(tree$edge[, 2], tree$edge[, 1]))

    linkD <- LinkDataFrame(nodeLab = fLab,
                       nodeLab_alias = faLab,
                       nodeNum = nd,
                       isLeaf = nd %in% leaf)


    # create the rowData
    return(linkD)
    }
