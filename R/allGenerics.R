
# -----------------------------------------------------------------------------
### Accessors for treeSummarizedExperiment
# -----------------------------------------------------------------------------
#' @importMethodsFrom SummarizedExperiment assays
#' @importFrom methods callNextMethod
#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("assays", signature("treeSummarizedExperiment"),
          function(x, use.nodeLab = FALSE, ..., withDimnames = FALSE){
              out <- callNextMethod(x, withDimnames)

              if (use.nodeLab) {
                  nodeLab <- x@linkData$nodeLab
                  if (any(duplicated(nodeLab))) {
                      nodeLab <- x@linkData$nodeLab_alias
                  }

                  outR <- lapply(out, function(x) {
                      rownames(x) <- nodeLab
                      x
                  })
              } else {
                  outR <- out
              }
              return(outR)
          })


#' @importMethodsFrom SummarizedExperiment rowData
#' @importFrom methods callNextMethod
#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("rowData", "treeSummarizedExperiment",
          function(x, use.names = TRUE, ...) {
              vv <- callNextMethod()
              cv <- unlist(lapply(vv, class))
              isInternal <- cv == "internal_rowData"
              vv[, !isInternal, drop = FALSE]
          })

#' @importMethodsFrom SummarizedExperiment "rowData<-"
#' @importFrom methods callNextMethod
#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("rowData", "treeSummarizedExperiment",
          function(x, ..., value) {
              callNextMethod()
          })

#' @rdname treeSummarizedExperiment-accessor
#' @export
setGeneric("linkData", function(x) {
    standardGeneric("linkData")
})


#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("linkData", signature("treeSummarizedExperiment"),
          function(x) {
              x@linkData
          })

#' @rdname treeSummarizedExperiment-accessor
#' @export
setGeneric("treeData", function(x) {
    standardGeneric("treeData")
})

#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("treeData", signature("treeSummarizedExperiment"),
          function(x) {
              x@treeData
          })

#' @importFrom methods callNextMethod
#' @importFrom SummarizedExperiment assays rowData colData
#  @importFrom BiocGenerics normalize
#' @importFrom S4Vectors metadata
#' @rdname treeSummarizedExperiment-accessor
#' @export
setMethod("[", signature(x = "treeSummarizedExperiment"),
          function(x, i, j){
              # Subset the traditional slots from SummarizedExperiment
              nx <- callNextMethod()

              # new slot
              linkD <- x@linkData
              if (!missing(i)) {
                  if (is.character(i)) {
                      fmt <- paste0("<", class(x),
                                    ">[i,] index out of bounds: %s")
                      i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                          i, rownames(x), fmt)}
                  i <- as.vector(i)
                  lk <- linkD[i, , drop = FALSE]
              }

              # update slots
              final <- BiocGenerics:::replaceSlots(nx,
                                                   linkData = lk)

              return(final)
          })


#' @keywords internal
#' @importFrom methods callNextMethod
#' @importMethodsFrom SummarizedExperiment show
setMethod("show", "treeSummarizedExperiment", function(object) {
    callNextMethod()
    cat(
        "treeData:", " a phylo \n",
        "linkData:", " a ", class(linkData(object)), " with ",
        ncol(linkData(object)), " columns \n",
        sep=""
    )
})



#' #' @keywords internal
#' setGeneric("show", function(x) {
#'     standardGeneric("show")
#' })


#' @keywords internal
#' @importMethodsFrom S4Vectors show
setMethod("show", "LinkDataFrame", function(object) {
    x_class <- class(object)
    left_len <- ncol(object@LinkData)
    right_len <- ncol(object)

    cat(x_class, " object with ",
        left_len, " link data ", ifelse(left_len == 1L, "column", "columns"),
        " and ",
        right_len, " metadata ", ifelse(right_len == 1L, "column", "columns"),
        ":\n", sep="")

    out <- S4Vectors:::makePrettyMatrixForCompactPrinting(object,
                                                          .nakedMatrix)
    classX <- c(lapply(object@LinkData, class), list("|" = " "), lapply(object, class))
    classX <- unlist(classX)

    classinfo <- S4Vectors:::makeClassinfoRowForCompactPrinting(object, classX)
    classinfo[1, "|"] <- ""

    ## A sanity check, but this should never happen!
    stopifnot(identical(colnames(classinfo), colnames(out)))

    out <- rbind(classinfo, out)

    print(out, quote=FALSE, right=TRUE, max=length(out))
})

#' @importFrom S4Vectors showAsCell
#' @keywords internal
.nakedMatrix <- function(x) {
    x_mcols <- x
    x_nmc <- ncol(x_mcols)
    ans <- as.matrix(x@LinkData)
    x_len <- max(nrow(x), nrow(ans))
    if (x_nmc > 0L) {
        tmp <- as.data.frame(lapply(x_mcols, showAsCell), optional=TRUE)
        ans <- cbind(ans, `|`=rep.int("|", x_len), as.matrix(tmp))

    } else {
        ans <-  cbind(ans, `|`=rep.int("|", x_len))

    }
    ans
}

# .show_LinkDataFrame <- function(x){
#
#     x_class <- class(x)
#     left_len <- ncol(x@LinkData)
#     right_len <- ncol(x)
#
#     cat(x_class, " object with ",
#         left_len, " link data ", ifelse(left_len == 1L, "column", "columns"),
#         " and ",
#         right_len, " metadata ", ifelse(right_len == 1L, "column", "columns"),
#         ":\n", sep="")
#
#     out <- S4Vectors:::makePrettyMatrixForCompactPrinting(x,
#                                                           .nakedMatrix)
#     classX <- c(lapply(x@LinkData, class), list("|" = " "), lapply(x, class))
#     classX <- unlist(classX)
#
#     classinfo <- S4Vectors:::makeClassinfoRowForCompactPrinting(x, classX)
#     classinfo[1, "|"] <- ""
#
#     ## A sanity check, but this should never happen!
#     stopifnot(identical(colnames(classinfo), colnames(out)))
#
#     out <- rbind(classinfo, out)
#
#     print(out, quote=FALSE, right=TRUE, max=length(out))
# }
