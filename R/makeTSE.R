#' A toy TreeSummarizedExperiment object
#' 
#' \code{makeTSE} creates a toy TreeSummarizedExperiment object. 
#' 
#' @param include.rowTree TRUE or FALSE. Default is TRUE, so the output
#'   \code{TreeSummarizedExperiment} has a \code{phylo} object in
#'   \code{rowTree}.
#' @param include.colTree TRUE or FALSE. Default is TRUE, so the output
#'   \code{TreeSummarizedExperiment} has a \code{phylo} object in
#'   \code{colTree}.
#' @param nrow a numeric value to specify the number of rows of
#'   \code{TreeSummarizedExperiment}
#' @param ncol a numeric value to specify the number of columns of
#'   \code{TreeSummarizedExperiment}
#' @importFrom ape rtree
#' @export
#' @return A TreeSummarizedExperiment object
#' @details The \code{assays} contains a matrix with values from
#'   \code{1:(nrow*ncol)}. The \code{rowData} has two columns, \code{var1} and
#'   \code{var2}. \code{var1} is created with \code{rep_len(letters, nrow)}.
#'   \code{var2} is created with \code{rep_len(c(TRUE, FALSE), nrow)}. The
#'   \code{colData} has two columns, \code{ID} and \code{group}. \code{ID} is
#'   created with \code{seq_len(ncol)}. \code{group} is created with
#'   \code{rep_len(LETTERS[1:2], ncol)}. The row/col tree is generated with
#'   \code{ape::rtree()}. So, to generate reproducible trees, \code{set.seed()}
#'   is required.
#' @author Ruizhu Huang
#' @examples 
#' 
#' set.seed(1)
#' makeTSE()

makeTSE <- function(nrow = 10, ncol = 4, 
                    include.rowTree = TRUE,
                    include.colTree = TRUE) {
    
    # the assays table
    tab <- matrix(seq_len(nrow*ncol), nrow = nrow)
    colnames(tab) <- paste0("sample", seq_len(ncol))
    rownames(tab) <- paste0("entity", seq_len(nrow))

    # the row data
    rD <- data.frame(var1 = rep_len(letters, nrow),
                     var2 = rep_len(c(TRUE, FALSE), nrow),
                     row.names = rownames(tab))
    # the column data
    cD <- data.frame(ID = seq_len(ncol),
                     group = rep_len(LETTERS[1:2], ncol),
                     row.names = colnames(tab))
    
    # the row & column tree
    rT <- rtree(nrow)
    rT$tip.label <- rownames(tab)
    cT <- rtree(ncol)
    cT$tip.label <- colnames(tab)
    
    if (!include.rowTree) {rT <- NULL}
    if (!include.colTree) {cT <- NULL}
    
    # conctruct TSE
    out <- TreeSummarizedExperiment(assays = list(tab),
                                    rowData = rD,
                                    colData = cD,
                                    rowTree = rT,
                                    colTree = cT)
    return(out)
}

