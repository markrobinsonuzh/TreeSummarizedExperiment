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
setGeneric("rowTree", function(x, whichTree = 1)
    standardGeneric("rowTree")
)


#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("rowTree", signature("TreeSummarizedExperiment"),
          function(x, whichTree = 1) {
              if (is.null(whichTree)) {return(x@rowTree)}
              x@rowTree[[whichTree]]
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("colTree", function(x, whichTree = 1)
    standardGeneric("colTree")
)


#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("colTree", signature("TreeSummarizedExperiment"),
          function(x, whichTree = 1) {
              if (is.null(whichTree)) {return(x@colTree)}
              x@colTree[[whichTree]]
          })


#' @importFrom methods callNextMethod
#' @importMethodsFrom SummarizedExperiment rbind
#' @rdname TreeSummarizedExperiment-combine
#' @export
#'
setMethod("rbind", signature = "TreeSummarizedExperiment",
          function(..., deparse.level = 1){
              # For slots inherited from SCE & SE
              nx <- callNextMethod()
              
              # a list of TSE
              args <- list(...)
              
              # Column tree & link data should be consistent in TSEs
              drop.colLinks  <- drop.rowLinks  <- FALSE
              isEq <- .is_equal_link(args, dim = "col")
              if (!isEq) {
                  warning("colTree & colLinks differ in the provided TSEs.",
                          "\n colTree & colLinks are dropped after 'rbind'")
                  drop.colLinks <- TRUE
              }
              
              # Row tree & link data should be all NULL or all non-NULL
              tList <- lapply(args, rowTree, whichTree = NULL)
              if (.any_null_in_list(tList) & !.all_null_in_list(tList)) {
                  warning("rowTree should be all NULL or non-NULL in TSEs.",
                          "\n rowTree & rowLinks are dropped after 'rbind'")
                  drop.rowLinks <- TRUE
              }
              
              final <- .combine_link_tree(x = nx, args = args, 
                                          drop.rowLinks = drop.rowLinks,
                                          drop.colLinks = drop.colLinks,
                                          bind = "rbind")
              return(final)
           })

#' @importFrom methods callNextMethod
#' @importMethodsFrom SummarizedExperiment rbind
#' @rdname TreeSummarizedExperiment-combine
#' @export
#'
setMethod("cbind", signature = "TreeSummarizedExperiment",
          function(..., deparse.level = 1){
              
              # For slots inherited from SCE & SE
              nx <- callNextMethod()
              
              # a list of TSE
              args <- list(...)
              
              # Row tree & link data should be consistent in TSEs
              drop.colLinks  <- drop.rowLinks  <- FALSE
              isEq <- .is_equal_link(args, dim = "row")
              if (!isEq) {
                  warning("rowTree & rowLinks differ in the provided TSEs.",
                          "\n rowTree & rowLinks are dropped after 'cbind'")
                  drop.rowLinks <- TRUE
              }
              
              # Row tree & link data should be all NULL or all non-NULL
              tList <- lapply(args, colTree, whichTree = NULL)
              if (.any_null_in_list(tList) & !.all_null_in_list(tList)) {
                  warning("colTree should be all NULL or non-NULL in TSEs.",
                          "\n colTree & colLinks are dropped after 'cbind'")
                  drop.colLinks <- TRUE
              }
              
              final <- .combine_link_tree(x = nx, args = args, 
                                          drop.rowLinks = drop.rowLinks,
                                          drop.colLinks = drop.colLinks,
                                          bind = "cbind")
              
              return(final)
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("referenceSeq", signature = c("x"),
           function(x) standardGeneric("referenceSeq"))

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("referenceSeq", signature = c(x = "TreeSummarizedExperiment"),
    function(x){
        x@referenceSeq
    }
)

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("referenceSeq<-", signature = c("x"),
           function(x, value) standardGeneric("referenceSeq<-"))

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setReplaceMethod("referenceSeq", signature = c(x = "TreeSummarizedExperiment"),
    function(x, value){
        x <- updateObject(x)
        if(is.null(value)){
         value <- value
        } else if(!is(value,"DNAStringSet") &&
               !is.list(value) &&
               !is(value,"DNAStringSetList")){
         value <- as(value,"DNAStringSet")
        } else if(!is(value,"DNAStringSetList") &&
               is.list(value)){
         value <- DNAStringSetList(value)
        }
        x <- .set_referenceSeq(x, value)
        validObject(x)
        x
    }
)

.set_referenceSeq <- function(x, value){
    x@referenceSeq <- value
    x
}


#' @importFrom methods callNextMethod
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("[", signature(x = "TreeSummarizedExperiment"),
          function(x, i, j, ..., drop = TRUE){
              x <- updateObject(x)
              # Subset the rowLinks
              lr <- rowLinks(x)
              rt <- rowTree(x)
              if (!missing(i) & !is.null(rt)) {
                  # match with rownames
                  # multiple rows in assays might have the same name
                  if (is.character(i)) {
                      isRn <- all(i %in% rownames(x))
                      if (isRn) {
                          i <- which(rownames(x) %in% i)
                      } else {
                          stop(i, " can't be found in rownames")
                      }
                  }

                  nlr <- lr[i, , drop = FALSE]
              } else {
                  nlr <- lr
              }

              # Subset the colLinks
              lc <- colLinks(x)
              ct <- colTree(x)
              if (!missing(j) & !is.null(ct)) {
                  # match with colnames
                  # multiple columns in assays might have the same name
                  if (is.character(j)) {
                      isCn <- all(j %in% colnames(x))
                      if (isCn) {
                          j <- which(colnames(x) %in% j)
                      } else {
                          stop(j, " can't be found in colnames")
                      }
                  }
                  nlc <- lc[j, , drop = FALSE]
              } else {
                  nlc <- lc
              }

              #
              refSeq <- referenceSeq(x)
              if (!missing(i)) {
                  refSeq <- referenceSeq(x)
                  if(!is.null(refSeq)){
                      if(is(refSeq,"DNAStringSetList")){
                          ii <- rep(list(i),length(refSeq))
                          refSeq <- refSeq[ii]
                      } else {
                          refSeq <- refSeq[i]
                      }
                  }
              }

              # Subset the traditional slots from SummarizedExperiment
              nx <- callNextMethod()

              # update slots
              final <- BiocGenerics:::replaceSlots(nx,
                                                   rowLinks = nlr,
                                                   colLinks = nlc,
                                                   referenceSeq = refSeq)
              validObject(final)
              return(final)
          }
)

#' @importFrom methods callNextMethod
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @rdname TreeSummarizedExperiment-accessor
#' @export
setReplaceMethod("[",
    signature(x = "TreeSummarizedExperiment", "ANY", "ANY", 
              "TreeSummarizedExperiment"),
    function(x, i, j, ..., value){
        x <- updateObject(x)
        value <- updateObject(value)
        if (missing(i) && missing(j)) {
            # do callNextMethod because of objects potentially being updated
            return(callNextMethod())
        }


        # TODO rowLinks
        # TODO colLinks

        if (!missing(i)) {
            x_refSeq <- referenceSeq(x)
            value_refSeq <- referenceSeq(value)
            if((!is.null(x_refSeq) & is.null(value_refSeq)) ||
               is.null(x_refSeq) & !is.null(value_refSeq) ||
               !is(x_refSeq, class(value_refSeq))){
                stop("'x' and 'value' must have the same type of ",
                     "referenceSeq()", call. = FALSE)
            }
            if(!is.null(x_refSeq)){
                if(is(x_refSeq,"DNAStringSetList")){
                    if(length(referenceSeq(value)) != length(x_refSeq)){
                        stop("DNAStringSetList as 'referenceSeq' must have ",
                             "the same length to be merged.", call. = FALSE)
                    }
                    ii <- rep(list(i),length(x_refSeq))
                    x_refSeq[ii] <- value_refSeq
                } else {
                    x_refSeq[i] <- value_refSeq
                }
            }
        }

        x <- callNextMethod()
        x <- BiocGenerics:::replaceSlots(x,
                                         referenceSeq = x_refSeq,
                                         check = FALSE)
        validObject(x)
        x
    }
)




#' @importFrom methods callNextMethod
#' @rdname TreeSummarizedExperiment-accessor
#' @export
setReplaceMethod("rownames", signature(x = "TreeSummarizedExperiment"),
                 function(x, value){
                     x <- updateObject(x)
                     x <- callNextMethod()
                     if(!is.null(x@rowLinks)){
                         rownames(x@rowLinks) <- value
                     }
                     x
                 }
)

#' @importFrom methods callNextMethod
#' @rdname TreeSummarizedExperiment-accessor
#' @export
setReplaceMethod("colnames", signature(x = "TreeSummarizedExperiment"),
                 function(x, value){
                     x <- updateObject(x)
                     x <- callNextMethod()
                     if(!is.null(x@colLinks)){
                         rownames(x@colLinks) <- value
                     }
                     x
                 }
)


#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("subsetByNode", function(x, rowNode, colNode)
    standardGeneric("subsetByNode")
)

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("subsetByNode", signature(x = "TreeSummarizedExperiment"),
          function(x, rowNode, colNode){
              x <- updateObject(x)
              # row link
              rl <- rowLinks(x)
              if (!missing(rowNode)) {
                  if (!is.numeric(rowNode)) {
                      rowNode <- convertNode(tree = rowTree(x), node = rowNode)
                  }
                  x <- x[which(rl$nodeNum %in% rowNode),]
              }

              # column link
              cl <- colLinks(x)
              if (!missing(colNode)) {
                  if (!is.numeric(colNode)) {
                      colNode <- convertNode(tree = colTree(x), node = colNode)
                  }
                  x <- x[, which(cl$nodeNum %in% colNode)]
              }
              return(x)
          }
)

#' @keywords internal
#' @importFrom methods callNextMethod
#' @importMethodsFrom SingleCellExperiment show
setMethod("show", "TreeSummarizedExperiment", function(object) {
    object <- updateObject(object)
    callNextMethod()

    rt <- rowTree(object, whichTree = NULL)
    ct <- colTree(object, whichTree = NULL)



    rlk <- rowLinks(object)
    clk <- colLinks(object)

    # on row
    if (is.null(rt)) {
        msg1a <- "rowLinks: NULL"
        msg1b <- "rowTree: NULL"
    } else {
        msg1a <- sprintf("rowLinks: a LinkDataFrame (%d %s)",
                         nrow(rlk), "rows")


        # the number of leaf nodes & internal nodes
        nlr <- sum(unlist(lapply(rt, countLeaf)))
        msg1b <- sprintf("rowTree: %d %s tree(s) (%d leaves)",
                         length(rt), class(rt[[1]]), nlr)
    }

    # on column
    if (is.null(ct)) {
        msg2a <- "colLinks: NULL"
        msg2b <- "colTree: NULL"
    } else {
        msg2a <- sprintf("colLinks: a LinkDataFrame (%d %s)", nrow(clk), "rows")

        # the number of leaf nodes & internal nodes
        nlc <- sum(unlist(lapply(ct, countLeaf)))
        msg2b <- sprintf("colTree: %d %s tree(s) (%d leaves)",
                         length(ct), class(ct[[1]]), nlc)

    }

    cat(msg1a, "\n", msg1b, "\n",
        msg2a, "\n", msg2b, "\n",
        sep = "")

    referenceSeq <- object@referenceSeq
    if(!is.null(referenceSeq)){
        if(is(referenceSeq,"DNAStringSetList")){
            msg <- sprintf(paste0("referenceSeq: a ", class(referenceSeq),
                                  " (%s x %s sequences each)\n"),
                           length(referenceSeq),
                           unique(lengths(referenceSeq)))
        } else {
            msg <- sprintf(paste0("referenceSeq: a ", class(referenceSeq),
                                  " (%s sequences)\n"),
                           length(referenceSeq))
        }
        cat(msg)
    }
})
