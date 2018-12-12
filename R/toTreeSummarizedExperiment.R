#' #' Translate a phyloseq object to a TreeSummarizedExperiment
#' #'
#' #' @param data A phyloseq object.
#' #'
#' #' @importFrom phyloseq tax_table phy_tree otu_table sample_data
#' #' @return a TreeSummarizedExperiment object
#' #' @author Ruizhu HUANG
#' #' @export
#' #' @examples
#' #'
#' #' taxTab <- data.frame(R1 = rep("A", 5),
#' #' R2 = c("B1", rep("B2", 4)),
#' #' R3 = c("C1", "C2", "C3", "C3", "C4"))
#' #'
#' #'
#' #' tree <- toTree(data = taxTab)
#' #'
#'
#' toTreeSummarizedExperiment <- function(data) {
#'     if (inherits(data, "phyloseq")) {
#'         # taxonomy
#'         tax <- tax_table(data)@.Data
#'
#'         # tree
#'         treeD <- phy_tree(data)
#'
#'         # count
#'         count <- otu_table(data)@.Data
#'
#'         # sample data
#'         sampD <- data.frame(sample_data(data))
#'         sampD <- DataFrame(sampD)
#'
#'         # output TreeSummarizedExperiment object
#'         tse <- treeSummarizedExperiment(tree = treeD,
#'                                         assays = list(count),
#'                                         rowData = tax,
#'                                         colData = sampD)
#'         return(tse)
#'     }
#' }
