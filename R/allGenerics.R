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
setGeneric("rowTree", function(x, whichTree = 1, value)
    standardGeneric("rowTree")
)

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("rowTree", signature("TreeSummarizedExperiment"),
          function(x, whichTree = 1, value) {
              if (is.null(whichTree)) {return(x@rowTree)}
              xx <- x@rowTree[whichTree]
              if (length(xx) == 1) {
                  xx <- xx[[1]]
              }
              return(xx)
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("rowTree<-", function(x, whichTree = 1, value)
    standardGeneric("rowTree<-")
)
#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("rowTree<-", signature("TreeSummarizedExperiment"),
          function(x, whichTree = 1, value) {
              if (is.null(value)) {
                  out <- BiocGenerics:::replaceSlots(x, 
                                                     rowTree = NULL, 
                                                     rowLinks = NULL)
                  return(out)
              }
              
              # 1) replace specified trees (e.g., whichTree = 1) 
              # 2) replace all trees (whichTree = NULL)
              
              out <- .replace_tree(x = x, value = value,
                                   whichTree = whichTree, 
                                   dim = "row")
              
              # new data
              drop <- out$drop
              if (length(drop)) {x <- x[-drop, ]}
              nlk <- out$new_links
              rownames(nlk) <- rownames(x)
              
              # update the row tree & link
              BiocGenerics:::replaceSlots(x, 
                                          rowTree = out$new_tree, 
                                          rowLinks = out$new_links)
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
              xx <- x@colTree[whichTree]
              if (length(xx) == 1) {
                  xx <- xx[[1]]
              }
              return(xx)
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("colTree<-", function(x, whichTree = 1, value)
    standardGeneric("colTree<-")
)
#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("colTree<-", signature("TreeSummarizedExperiment"),
          function(x, whichTree = 1, value) {
              
              if (is.null(value)) {
                  out <- BiocGenerics:::replaceSlots(x, 
                                                     colTree = NULL, 
                                                     colLinks = NULL)
                  return(out)
              }
              
              # 1) replace specified trees (e.g., whichTree = 1) 
              # 2) replace all trees (whichTree = NULL)
              
              out <- .replace_tree(x = x, value = value,
                                   whichTree = whichTree, 
                                   dim = "col")
              # new data
              drop <- out$drop
              if (length(drop)) {x <- x[, -drop]}
              nlk <- out$new_links
              rownames(nlk) <- colnames(x)
              
              # update the col tree & link
              BiocGenerics:::replaceSlots(x, 
                                          colTree = out$new_tree, 
                                          colLinks = out$new_links)
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("rowTreeNames", function(x, value)
    standardGeneric("rowTreeNames")
)

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("rowTreeNames", signature("TreeSummarizedExperiment"),
          function(x, value) {
              rT <- rowTree(x, whichTree = NULL)
              names(rT)
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("rowTreeNames<-", signature = c("x"),
           function(x, value) standardGeneric("rowTreeNames<-"))

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("rowTreeNames<-", signature("TreeSummarizedExperiment"),
          function(x, value) {
              # the row tree
              rT <- rowTree(x, whichTree = NULL)
              namePair <- setNames(value, names(rT))
              names(rT) <- value
              
              # the row link
              rL <- rowLinks(x)
              
              
              # update the column whichTree in the rowLinks
              nrL <- .update_whichTree(rL, namePair)
              
              # update the row tree & link
              BiocGenerics:::replaceSlots(x, rowTree = rT, rowLinks = nrL)
          })



#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("colTreeNames", function(x, value)
    standardGeneric("colTreeNames")
)

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("colTreeNames", signature("TreeSummarizedExperiment"),
          function(x, value) {
              cT <- colTree(x, whichTree = NULL)
              names(cT)
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("colTreeNames<-", 
           function(x, value) standardGeneric("colTreeNames<-"))

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("colTreeNames<-", signature("TreeSummarizedExperiment"),
          function(x, value) {
              # the column tree
              cT <- colTree(x, whichTree = NULL)
              namePair <- setNames(value, names(cT))
              names(cT) <- value
              
              # the column link
              cL <- colLinks(x)
              
              
              # update the column whichTree in the colLinks
              ncL <- .update_whichTree(cL, namePair)
              
              # update the row tree & link
              BiocGenerics:::replaceSlots(x, colTree = cT, colLinks = ncL)
          })



#' @importFrom methods callNextMethod
#' @importMethodsFrom SummarizedExperiment rbind
#' @rdname TreeSummarizedExperiment-combine
#' @export
setMethod("rbind", signature = "TreeSummarizedExperiment",
          function(..., deparse.level = 1){
              
              old <- S4Vectors:::disableValidity()
              if (!isTRUE(old)) {
                  S4Vectors:::disableValidity(TRUE)
                  on.exit(S4Vectors:::disableValidity(old))
              }
              
              # For slots inherited from SCE & SE
              nx <- callNextMethod()
              
              # For slots only in TSE:
              # ------------------------------------------------------------
              args <- list(...)
              
              # Column tree & link data should be consistent in TSEs
              drop.colLinks  <- drop.rowLinks  <- FALSE
              isEq <- .is_equal_link(args, dim = "col")
              if (!isEq) {
                  warning("colTree & colLinks differ in the provided TSEs.",
                          "\n colTree & colLinks are dropped after 'rbind'", 
                          call. = FALSE)
                  drop.colLinks <- TRUE
              }
              
              # Row tree & link data should be all NULL or all non-NULL
              tList <- lapply(args, rowTree, whichTree = NULL)
              if (.any_null_in_list(tList) & !.all_null_in_list(tList)) {
                  warning("rowTree should be all NULL or non-NULL in TSEs.",
                          "\n rowTree & rowLinks are dropped after 'rbind'",
                          call. = FALSE)
                  drop.rowLinks <- TRUE
              }
              
              nnx <- .bind_link_tree(x = nx, args = args,
                                     drop.rowLinks = drop.rowLinks,
                                     drop.colLinks = drop.colLinks,
                                     bind = "rbind")
              
              # rbind on the referenceSeq slot
              refSeq <- .rbind_refSeq(args)
              BiocGenerics:::replaceSlots(nnx,
                                          referenceSeq = refSeq,
                                          check = FALSE)
          })

#' @importFrom methods callNextMethod
#' @importMethodsFrom SummarizedExperiment rbind
#' @rdname TreeSummarizedExperiment-combine
#' @export
#'
setMethod("cbind", signature = "TreeSummarizedExperiment",
          function(..., deparse.level = 1){
              
              old <- S4Vectors:::disableValidity()
              if (!isTRUE(old)) {
                  S4Vectors:::disableValidity(TRUE)
                  on.exit(S4Vectors:::disableValidity(old))
              }
              
              # For slots inherited from SCE & SE
              nx <- callNextMethod()
              
              # For slots only in TSE:
              # ------------------------------------------------------------
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
              
              .bind_link_tree(x = nx, args = args,
                                 drop.rowLinks = drop.rowLinks,
                                 drop.colLinks = drop.colLinks,
                                 bind = "cbind")
              
              
          })

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setGeneric("referenceSeq", signature = c("x"),
           function(x) standardGeneric("referenceSeq"))

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("referenceSeq", signature = c(x = "TreeSummarizedExperiment"),
    function(x){
        seq <- x@referenceSeq
        if(!is.null(seq)){
            if(is(seq,"DNAStringSetList")){
                seq_u <- unlist(seq)
                names(seq_u) <- unlist(rep(rownames(x),length(seq)))
                seq <- relist(seq_u,seq)
            } else {
                names(seq) <- rownames(x)
            }
        }
        seq
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
    if(is(value,"DNAStringSetList")){
        value <- relist(unname(unlist(value)),value)
    } else {
        value <- unname(value)
    }
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
              
              # Subset rowLinks & rowTree & referenceSeq
              nlr <- lr <- rowLinks(x)
              nrt <- rt <- rowTree(x, whichTree = NULL)
              refSeq <- referenceSeq(x)
              
              if (!missing(i)) {
                  i <- .numeric_ij(ij = i, x = x, dim = "row")
                  
                  if(!length(i) | !sum(i)) {
                      refSeq <- nlr <- nrt <- NULL
                  } else {
                      # referenceSeq
                      if(!is.null(refSeq)){
                          if(is(refSeq,"DNAStringSetList")){
                              if (any(i < 0)) {
                                  abs_i <- abs(i)
                                  li <- lapply(lengths(refSeq), seq_len)
                                  ii <- lapply(li, FUN = function(x) {
                                      !x %in% abs_i
                                  })
                              } else {
                                  ii <- rep(list(i),length(refSeq))  
                              }
                              refSeq <- refSeq[ii]
                          } else {refSeq <- refSeq[i]}
                      }
                      # rowLinks & rowTree
                      if (!is.null(rt)) {
                          nlr <- lr[i, , drop = FALSE]
                          nrt <- rt[unique(nlr$whichTree)]
                      }
                  }
              }
              

              # Subset the colLinks & colTree
              nlc <- lc <- colLinks(x)
              nct <- ct <- colTree(x, whichTree = NULL)
              if (!missing(j)) {
                  j <- .numeric_ij(ij = j, x = x, dim = "col")
                  if(!length(j) | !sum(j)) {
                      nlc <- nct <- NULL
                  } else {
                      if (!is.null(ct)) {
                          nlc <- lc[j, , drop = FALSE]
                          nct <- ct[unique(nlc$whichTree)]
                      }
                  }
              }

              # Subset the traditional slots from SummarizedExperiment
              nx <- callNextMethod()

              # update slots
              final <- BiocGenerics:::replaceSlots(nx,
                                                   rowLinks = nlr,
                                                   colLinks = nlc,
                                                   rowTree = nrt,
                                                   colTree = nct,
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
                     
                     # TODO rowLinks and rowTree
                     # TODO colLinks
                     if ((sum(missing(i), missing(j)) == 1)) {
                         outR <- .replace_link_tree_1d(x = x, value = value, ij = i,
                                                       dim = "row")
                         outC <- .replace_link_tree_1d(x = x, value = value,
                                                       ij = j, dim = "col")
                     } else {
                         out <- .replace_link_tree_2d(x = x, value = value, 
                                                      i = i, j = j)
                         outR <- out$outR
                         outC <- out$outC
                     }
                     
                     
                     x_refSeq <- referenceSeq(x)
                     if (!missing(i)) {
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
                                 if (any(i < 0)) {
                                     abs_i <- abs(i)
                                     li <- lapply(lengths(x_refSeq), seq_len)
                                     ii <- lapply(li, FUN = function(x) {
                                         !x %in% abs_i
                                     })
                                 } else {
                                     ii <- rep(list(i),length(x_refSeq))  
                                 }
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
                     if (!is.null(outR)) {
                         nl <- outR$new_links
                         rownames(nl) <- rownames(x)
                         x <- BiocGenerics:::replaceSlots(x,
                                                          rowTree = outR$new_tree,
                                                          rowLinks = nl,
                                                          check = FALSE)
                     }
                     
                     if (!is.null(outC)) {
                         nl <- outC$new_links
                         rownames(nl) <- colnames(x)
                         x <- BiocGenerics:::replaceSlots(x,
                                                          colTree = outC$new_tree,
                                                          colLinks = nl,
                                                          check = FALSE)
                     }
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
setGeneric("subsetByNode", function(x, rowNode, colNode,
                                    whichRowTree, whichColTree)
    standardGeneric("subsetByNode")
)

#' @rdname TreeSummarizedExperiment-accessor
#' @export
setMethod("subsetByNode", signature(x = "TreeSummarizedExperiment"),
          function(x, rowNode, colNode, whichRowTree, whichColTree){
              x <- updateObject(x)
              # row link
              rl <- rowLinks(x)
              rt <- rowTree(x, whichTree = NULL)
              if (!missing(whichRowTree)) {
                  rt <- rt[whichRowTree]
                  rnam <- names(rt)
                  x <- x[rl$whichTree %in% rnam,]
                  x <- BiocGenerics:::replaceSlots(object = x, 
                                                   rowTree = rt)
                  rl <- rowLinks(x)
              }
              
              if (!missing(whichRowTree) |
                  !missing(rowNode)) {
                  if (!length(rt)) {
                      warning("The row tree is not available.", 
                              call. = FALSE)}
              }
              
              if (!missing(rowNode)) {
                  if (!is.numeric(rowNode)) {
                      rowNode <- unlist(lapply(rt, convertNode, node = rowNode))
                  }
                  x <- x[which(rl$nodeNum %in% rowNode),]
              }
              
              
              
              # column link
              cl <- colLinks(x)
              ct <- colTree(x, whichTree = NULL)
              if (!missing(whichColTree)) {
                  ct <- ct[whichColTree]
                  cnam <- names(ct)
                  x <- x[, cl$whichTree %in% cnam]
                  x <- BiocGenerics:::replaceSlots(object = x, 
                                                   colTree = ct)
                  cl <- colLinks(x)
              }
              
              if (!missing(whichColTree) |
                  !missing(colNode)) {
                  if (!length(ct)) {
                      warning("The column tree is not available.", 
                              call. = FALSE)}
              }
              
              if (!missing(colNode)) {
                  if (!is.numeric(colNode)) {
                      colNode <- unlist(lapply(ct, convertNode, node = colNode))
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
