.sce_to_tse <- function(sce){
    tse <- new("TreeSummarizedExperiment",
               sce,
               rowTree = NULL,
               colTree = NULL,
               rowLinks = NULL,
               colLinks = NULL)
    tse
}

#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setAs("SingleCellExperiment", "TreeSummarizedExperiment",function(from) {
    .sce_to_tse(from)
})

#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setAs("SummarizedExperiment", "TreeSummarizedExperiment",function(from) {
    .sce_to_tse(as(from,"SingleCellExperiment"))
})

#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
setAs("RangedSummarizedExperiment", "TreeSummarizedExperiment",function(from) {
    .sce_to_tse(as(from,"SingleCellExperiment"))
})