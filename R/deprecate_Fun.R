#' @rdname findDescendant
#' @export
findOS <- function(tree,
                   node,
                   only.leaf = TRUE,
                   self.include = FALSE,
                   use.alias = FALSE){
    .Deprecated("findDescendant")
    findDescendant(tree = tree,
                   node = node,
                   only.leaf = only.leaf,
                   self.include = self.include,
                   use.alias = use.alias)
}

#' @rdname convertNode
#' @export
transNode <- function(tree, node, use.alias = FALSE,
                      message = FALSE){
    .Deprecated("convertNode")
    convertNode(tree = tree, node = node, use.alias = use.alias,
                message = message)
}

#' @rdname aggTSE
#' @export
aggValue <- function(x, rowLevel = NULL, rowBlock = NULL, 
                     colLevel = NULL, colBlock = NULL, 
                     FUN = sum, assay = NULL,
                     message = FALSE){
    .Deprecated("aggTSE")
    aggTSE(x = x, 
           rowLevel = rowLevel, rowBlock = rowBlock, 
           colLevel = colLevel, colBlock = colBlock, 
           FUN = sum, whichAssay = assay,
           message = message,
           whichRowTree = 1, whichColTree = 1,
           rowFirst = TRUE)
}