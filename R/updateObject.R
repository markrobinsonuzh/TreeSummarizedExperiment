#' Update a \code{TreeSummarizedExperiment} object
#' 
#' Update \code{TreeSummarizedExperiment} objects to the latest version of the 
#' class structure. This is usually called by methods in the 
#' \code{TreeSummarizedExperiment} package rather than by users or downstream 
#' packages.
#' 
#' @param object A \code{TreeSummarizedExperiment} object
#' 
#' @param ... additional arguments, for use in specific \code{updateObject}
#'   methods.
#'   
#' @param verbose \code{TRUE} or \code{FALSE}, indicating whether information 
#'   about the update should be reported.
#'   
#' @export
#'   
#' @return An updated \code{TreeSummarizedExperiment} object
setMethod("updateObject", "TreeSummarizedExperiment",
    function(object, ..., verbose = FALSE){
        old <- S4Vectors:::disableValidity()
        if (!isTRUE(old)) {
            S4Vectors:::disableValidity(TRUE)
            on.exit(S4Vectors:::disableValidity(old))
        }
        object <- callNextMethod()
        if(!.hasSlot(object,"referenceSeq")){
            if(verbose){
                message("Updating ", class(object)[1], " object ...\n",
                        appendLF = FALSE)
            }
            object@referenceSeq <- NULL
        }
        object
    }
)
