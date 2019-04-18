# -----------------------------------------------------------------------------
### Accessors for TreeSummarizedExperiment
# -----------------------------------------------------------------------------

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("rowLink", function(x)
    standardGeneric("rowLink")
)


#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("rowLink", signature("TreeSummarizedExperiment"),
          function(x) {
              x@rowLink
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("colLink", function(x)
    standardGeneric("colLink")
)


#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("colLink", signature("TreeSummarizedExperiment"),
          function(x) {
              x@colLink
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

#' @keywords internal
#' @importFrom methods callNextMethod
#' @importMethodsFrom SingleCellExperiment show
setMethod("show", "TreeSummarizedExperiment", function(object) {
    callNextMethod()

    rt <- rowTree(object)
    ct <- colTree(object)



    rlk <- rowLink(object)
    clk <- colLink(object)

    # on row
    if (is.null(rt)) {
        msg1a <- "rowLink:"
        msg1b <- "rowTree:"
    } else {
        msg1a <- sprintf("rowLink: a LinkDataFrame (%d %s)",
                         nrow(rlk), "rows")


        # the number of leaf nodes & internal nodes
        nlr <- countLeaf(rt)
        nnr <- countNode(rt) - countLeaf(rt)
        msg1b <- sprintf("rowTree: a %s (%d %s, %d %s)", class(rt),
                         nlr, "leaves", nnr, "internal nodes")
    }

    # on column
    if (is.null(ct)) {
        msg2a <- "colLink:"
        msg2b <- "colTree:"
    } else {
        msg2a <- sprintf("colLink: a LinkDataFrame (%d %s)", nrow(clk), "rows")

        # the number of leaf nodes & internal nodes
        nlc <- countLeaf(ct)
        nnc <- countNode(ct) - countLeaf(ct)
        msg2b <- sprintf("colTree: a %s (%d %s, %d %s)", class(ct),
                         nlc, "leaves", nnc, "internal nodes")
    }

    cat(msg1a, "\n", msg1b, "\n",
        msg2a, "\n", msg2b, "\n",
        sep = "")
})








