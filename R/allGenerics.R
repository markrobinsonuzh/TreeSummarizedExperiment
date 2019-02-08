
# -----------------------------------------------------------------------------
### Accessors for TreeSummarizedExperiment
# -----------------------------------------------------------------------------

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("linkData", function(x, onRow = TRUE) {
    standardGeneric("linkData")
})


#' @rdname TreeSummarizedExperiment-accessor
#' @importFrom methods is
#' @export
setMethod("linkData", signature("TreeSummarizedExperiment"),
          function(x, onRow = TRUE) {
              if (onRow) {
                  rD <- rowData(x)

                  # extract LinkData if LinkDataFrame, otherwise return NULL
                  if (is(rD, "LinkDataFrame")) {
                      rD@LinkData
                  } else { NULL }

              } else {
                  cD <- colData(x)

                  # extract LinkData if LinkDataFrame, otherwise return NULL
                  if (is(cD, "LinkDataFrame")) {
                      cD@LinkData
                  } else { NULL }
              }
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("treeData", function(x, onRow = TRUE) {
    standardGeneric("treeData")
})

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("treeData", signature("TreeSummarizedExperiment"),
          function(x, onRow = TRUE) {
              if (onRow) {
                  tD <- x@treeData$rowTree

              } else {
                  tD <- x@treeData$colTree
              }

              tD
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("metaCol", function(x, onRow = TRUE) {
    standardGeneric("metaCol")
})

#' @rdname TreeSummarizedExperiment-accessor
#' @importFrom methods is
#' @importFrom S4Vectors DataFrame
#' @export
setMethod("metaCol", signature("TreeSummarizedExperiment"),
          function(x, onRow = TRUE) {
              if (onRow) {

                  rD <- rowData(x)
                  if (is(rD, "LinkDataFrame")) {
                      rD1 <- as.data.frame(rD)
                      rD2 <- rD1[, colnames(rD1) %in% colnames(rD)]
                      DataFrame(rD2)
                  } else { rD }

              } else {

                  cD <- colData(x)
                  if (is(cD, "LinkDataFrame")) {
                      cD1 <- as.data.frame(cD)
                      cD2 <- cD1[, colnames(cD1) %in% colnames(cD)]
                      DataFrame(cD2)
                  } else { cD }

              }})
#' @keywords internal
#' @importFrom methods callNextMethod
#' @importMethodsFrom SummarizedExperiment show
setMethod("show", "TreeSummarizedExperiment", function(object) {
    callNextMethod()

    rt <- treeData(object, onRow = TRUE)
    ct <- treeData(object, onRow = FALSE)


    if (is.null(rt)) {
        msg1 <- ""
        nt1 <- 0
    } else {
        msg1 <- "exists in rowData. \n"
        nt1 <- 1
        }

    if (is.null(ct)) {
        msg2 <- ""
        nt2 <- 0
    } else {
        msg2 <- "exists in colData. \n"
        nt2 <- 1
        }

    cat(
        "treeData: ", nt1, " row tree; ", nt2, " column tree. \n",
        "linkData: ", msg1, msg2,
        sep=""
    )
})


# -----------------------------------------------------------------------------
### Accessors for LinkDataFrame
# -----------------------------------------------------------------------------
#' @importFrom methods callNextMethod
#' @rdname LinkDataFrame-accessor
#' @export
setMethod("[", signature(x = "LinkDataFrame"),
          function(x, i, j){

              # Subset the slots from DataFrame
              nx <- callNextMethod(x, i, j, drop = FALSE)

              # subset the new slot LinkData: only on rows
              nc <- ncol(x@LinkData)
              jj <- seq_len(nc)

              y <- x@LinkData
              ld <- callNextMethod(y, i, jj)

              # update slots
              final <- BiocGenerics:::replaceSlots(nx,
                                                   LinkData = ld)

              return(final)
              })

#' @importFrom methods callNextMethod
#' @rdname LinkDataFrame-accessor
#' @export
setMethod("as.data.frame", "LinkDataFrame", function(x) {
    nx <- callNextMethod(x)
    y <- x@LinkData
    nl <- callNextMethod(y)

    if (sum(dim(nx))) {
        final <- cbind(nl, nx)
    } else {
        final <- nl
    }

    return(final)
})


#' @importFrom methods callNextMethod
#' @rdname LinkDataFrame-accessor
#' @export
setMethod("$", signature(x = "LinkDataFrame"),
          function(x, name){
              nx <- callNextMethod()
              return(nx)
 })

#' @importFrom methods callNextMethod
#' @rdname LinkDataFrame-accessor
#' @export
setMethod("$<-", signature(x = "LinkDataFrame"),
          function(x, name, value){
              nx <- callNextMethod()
              # update slots
              final <- BiocGenerics:::replaceSlots(nx,
                                                   LinkData = x@LinkData)
              return(final)
          })

#' The show method of the class \strong{LinkDataFrame} is modified based on the
#' codes of the show method of the class \strong{"GPos"}. The original code
#' could be found in the
#' https://github.com/Bioconductor/GenomicRanges/blob/master/R/GPos-class.R
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

