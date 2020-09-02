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