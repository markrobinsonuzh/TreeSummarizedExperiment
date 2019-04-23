# -----------------------------------------------------------------------------
### Accessors for TreeSummarizedExperiment
# -----------------------------------------------------------------------------

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("rowLinks", function(x)
    standardGeneric("rowLinks")
)


#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("rowLinks", signature("TreeSummarizedExperiment"),
          function(x) {
              x@rowLinks
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("colLinks", function(x)
    standardGeneric("colLinks")
)


#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("colLinks", signature("TreeSummarizedExperiment"),
          function(x) {
              x@colLinks
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("rowTree", function(x)
    standardGeneric("rowTree")
)


#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("rowTree", signature("TreeSummarizedExperiment"),
          function(x) {
              x@rowTree
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("colTree", function(x)
    standardGeneric("colTree")
)


#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("colTree", signature("TreeSummarizedExperiment"),
          function(x) {
              x@colTree
          })


#' @importFrom methods callNextMethod
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @rdname TreeSummarizedExperiment-accessor
#' @export
#'
setMethod("[", signature(x = "TreeSummarizedExperiment"),
          function(x, i, j, ..., drop = TRUE){

             # Subset the traditional slots from SummarizedExperiment
              nx <- callNextMethod()

              # Subset the rowLinks
              lr <- x@rowLinks
              rt <- x@rowTree
              if (!missing(i) & !is.null(rt)) {
                  nlr <- lr[i, , drop = FALSE]
              } else {
                  nlr <- lr
              }

              # Subset the colLinks
              lc <- x@colLinks
              ct <- x@colTree
              if (!missing(j) & !is.null(ct)) {
                  nlc <- lc[j, , drop = FALSE]
              } else {
                  nlc <- lc
              }


              # update slots
              final <- BiocGenerics:::replaceSlots(nx,
                                                   rowLinks = nlr,
                                                   colLinks = nlc)

              return(final)
          })

#' @keywords internal
#' @importFrom methods callNextMethod
#' @importMethodsFrom SingleCellExperiment show
setMethod("show", "TreeSummarizedExperiment", function(object) {
    callNextMethod()

    rt <- rowTree(object)
    ct <- colTree(object)



    rlk <- rowLinks(object)
    clk <- colLinks(object)

    # on row
    if (is.null(rt)) {
        msg1a <- "rowLinks:"
        msg1b <- "rowTree:"
    } else {
        msg1a <- sprintf("rowLinks: a LinkDataFrame (%d %s)",
                         nrow(rlk), "rows")


        # the number of leaf nodes & internal nodes
        nlr <- countLeaf(rt)
        nnr <- countNode(rt) - countLeaf(rt)
        msg1b <- sprintf("rowTree: a %s (%d leaves)", class(rt), nlr)
    }

    # on column
    if (is.null(ct)) {
        msg2a <- "colLinks:"
        msg2b <- "colTree:"
    } else {
        msg2a <- sprintf("colLinks: a LinkDataFrame (%d %s)", nrow(clk), "rows")

        # the number of leaf nodes & internal nodes
        nlc <- countLeaf(ct)
        nnc <- countNode(ct) - countLeaf(ct)
        msg2b <- sprintf("colTree: a %s (%d leaves)", class(ct), nlc)
    }

    cat(msg1a, "\n", msg1b, "\n",
        msg2a, "\n", msg2b, "\n",
        sep = "")
})








