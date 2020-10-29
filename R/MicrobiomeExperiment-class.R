#' @include allClass.R
NULL

#' The \code{MicrobiomeExperiment} class
#'
#' The \code{MicrobiomeExperiment} class is a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}-like
#' class for microbiome data inheriting from
#' \code{\link[=TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' and
#' \code{\link[SingleCellExperiment:SingleCellExperiment]{SingleCellExperiment}}.
#'
#' \code{MicrobiomeExperiment} adds to the above mentioned class a slot for
#' storing a reference sequence as a
#' \code{\link[Biostrings:XStringSet-class]{DNAStringSet}}.
#'
#' @param ... Arguments passed to
#'   \code{\link[=TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' @param referenceSeq a \code{DNAStringSet} object or some object
#'   coercible to a \code{DNAStringSet} object. See
#'   \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} for more
#'   details.
#'
#' @name MicrobiomeExperiment-class
#'
#' @importClassesFrom Biostrings XStringSet DNAStringSet XStringSetList
#'   DNAStringSetList
#'
#' @importFrom Biostrings DNAStringSet DNAStringSetList
#'
#' @examples
#' sampleNames <- letters[1:4]
#' cd <- DataFrame(a = letters[1:4], b = 1:4)
#' rd <- DataFrame(Kingdom = "Bacteria",
#'                 Phylum = c("Firmicutes","Firmicutes","Firmicutes",
#'                            "Bacteroidetes","Euryarchaeota"),
#'                 Class = c("Bacilli","Bacilli","Bacilli","Saprospirae",
#'                           "Methanobacteria"),
#'                 Order = c("Bacillales","Lactobacillales","Lactobacillales",
#'                           "Saprospirales","Methanobacteriales"),
#'                 Family = c("Planococcaceae","Enterococcaceae","Enterococcaceae",
#'                            "Chitinophagaceae","Methanobacteriaceae"),
#'                 Genus = c("Staphylococcus",NA,"Melissococcus",
#'                           "Sediminibacterium","Methanobrevibacter"))
#' counts <- matrix(sample(1:1000, nrow(rd) * nrow(cd), replace=TRUE),
#'                  nr = nrow(rd),
#'                  nc = nrow(cd))
#' refSeq <- DNAStringSetList(one = DNAStringSet(c("A","A","A","A","A")),
#'                            two = DNAStringSet(c("A","A","A","A","A")))
#' me <- MicrobiomeExperiment(assays = SimpleList(counts = counts),
#'                            rowData = rd,
#'                            colData = cd,
#'                            referenceSeq = refSeq)
#' me
NULL

setClassUnion("DNAStringSetList_OR_DNAStringSet_OR_NULL",
              c("DNAStringSetList", "DNAStringSet","NULL"))

#' @rdname MicrobiomeExperiment-class
#' @export
setClass("MicrobiomeExperiment",
         contains = "TreeSummarizedExperiment",
         slots = c(referenceSeq = "DNAStringSetList_OR_DNAStringSet_OR_NULL"),
         prototype = list(referenceSeq = NULL)
)

#' @rdname MicrobiomeExperiment-internal
#' @export
setMethod("vertical_slot_names", "MicrobiomeExperiment",
          function(x) c("referenceSeq", callNextMethod())
)

################################################################################
# validity

.valid.MicrobiomeExperiment <- function(x)
{
    x_nrow <- length(x)
    if(!is.null(x@referenceSeq)){
        if(is(x@referenceSeq, "DNAStringSet")){
            referenceSeq_len <- length(x@referenceSeq)
        } else if(is(x@referenceSeq, "DNAStringSetList")){
            referenceSeq_len <- lengths(x@referenceSeq)
            referenceSeq_len <- unique(referenceSeq_len)
        }
        if (length(referenceSeq_len) != 1L) {
            txt <- "\n  lengths of 'referenceSeq' must all be equal."
            return(txt)
        }
        if (referenceSeq_len != x_nrow) {
            txt <- sprintf(
                paste0("\n  length(s) of 'referenceSeq' (%d) must equal nb of ",
                       "rows in 'x' (%d)"),
                referenceSeq_len, x_nrow)
            return(txt)
        }
    }
    NULL
}

S4Vectors::setValidity2("MicrobiomeExperiment", .valid.MicrobiomeExperiment)

################################################################################
# constructor

#' @rdname MicrobiomeExperiment-class
#' @export
MicrobiomeExperiment <- function(..., referenceSeq = NULL) {
    tse <- TreeSummarizedExperiment(...)
    .tse_to_me(tse, referenceSeq)
}

.tse_to_me <- function(tse, referenceSeq = NULL){
    me <- new("MicrobiomeExperiment",
              tse,
              referenceSeq = referenceSeq)
    me
}

################################################################################
# coercion

setAs("TreeSummarizedExperiment", "MicrobiomeExperiment", function(from) {
    .tse_to_me(from)
})

#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
setAs("RangedSummarizedExperiment", "MicrobiomeExperiment", function(from) {
    .tse_to_me(as(from,"TreeSummarizedExperiment"))
})

#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setAs("SummarizedExperiment", "MicrobiomeExperiment", function(from) {
    .tse_to_me(as(from,"TreeSummarizedExperiment"))
})

#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setAs("SingleCellExperiment", "MicrobiomeExperiment", function(from) {
    .tse_to_me(as(from,"TreeSummarizedExperiment"))
})

################################################################################
# accessors

#' Microbiome data methods
#'
#' Methods to get or set reference sequence data on a
#' \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}} object.
#'
#' @param x a \code{\link[=MicrobiomeExperiment-class]{MicrobiomeExperiment}}
#'   object
#'
#' @param value a a \code{\link[Biostrings:XStringSet-class]{DNAStringSet}}
#'   object or an object coercible to one.
#'
#' @name referenceSeq
#'
#' @export
setGeneric("referenceSeq", signature = c("x"),
           function(x) standardGeneric("referenceSeq"))
#' @rdname referenceSeq
#' @export
setMethod("referenceSeq", signature = c(x = "MicrobiomeExperiment"),
    function(x){
        x@referenceSeq
    }
)

#' @rdname referenceSeq
#' @export
setGeneric("referenceSeq<-", signature = c("x"),
           function(x, value) standardGeneric("referenceSeq<-"))
#' @rdname referenceSeq
#' @export
setReplaceMethod("referenceSeq", signature = c(x = "MicrobiomeExperiment"),
    function(x, value){
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

################################################################################
# subsetting

#' @rdname MicrobiomeExperiment-internal
setMethod("[", signature = c("MicrobiomeExperiment", "ANY", "ANY"),
    function(x, i, j, ..., drop = TRUE) {
        if (!missing(i)) {
            ans_refSeq <- referenceSeq(x)
            if(!is.null(ans_refSeq)){
                if(is(ans_refSeq,"DNAStringSetList")){
                    ii <- rep(list(i),length(ans_refSeq))
                    ans_refSeq <- ans_refSeq[ii]
                } else {
                    ans_refSeq <- ans_refSeq[i]
                }
                x <- BiocGenerics:::replaceSlots(x, ...,
                                                 referenceSeq = ans_refSeq,
                                                 check = FALSE)
            }
        }

        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(FALSE))
        x <- callNextMethod()
        S4Vectors:::disableValidity(FALSE)
        validObject(x)
        x
    }
)

#' @rdname MicrobiomeExperiment-internal
setReplaceMethod("[", signature = c("MicrobiomeExperiment", "ANY", "ANY", "MicrobiomeExperiment"),
    function(x, i, j, ..., value) {
        if (missing(i) && missing(j)) {
            return(value)
        }

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
                x <- BiocGenerics:::replaceSlots(x, ...,
                                                 referenceSeq = x_refSeq,
                                                 check = FALSE)
            }
        }

        S4Vectors:::disableValidity(TRUE)
        on.exit(S4Vectors:::disableValidity(FALSE))
        x <- callNextMethod()
        S4Vectors:::disableValidity(FALSE)
        validObject(x)
        x
    }
)


################################################################################
# show

setMethod("show", signature = c(object = "MicrobiomeExperiment"),
    function(object){
        callNextMethod()
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
        } else {
            msg <- "referenceSeq: NULL\n"
        }
        cat(msg)
    }
)
