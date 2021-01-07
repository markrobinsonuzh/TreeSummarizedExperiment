#' test all elements in a list are equal
#' @keywords internal
#' @examples 
#' \dontrun{
#' l1 <- list(a = 1, b = 2, c = 3)
#' l2 <- list(a = 1, b = 1, c = 1)
#' .all_equal_in_list(l1)
#' .all_equal_in_list(l2)
#' }
.all_equal_in_list <- function(x) {
    ux <- unique(x)
    length(ux) == 1
}

#' all elements in the list are NULL
#' @keywords internal
.all_null_in_list <- function(x) {
    xl <- lapply(x, is.null)
    all(unlist(xl))
}

#' all elements in the list are NULL
#' @keywords internal
.all_nonnull_in_list <- function(x) {
    xl <- lapply(x, !is.null)
    all(unlist(xl))
}

#' name y with x
#' @keywords internal
#' @examples 
#' \dontrun{
#' x <- letters[1:5]
#' y <- 1:5
#' .name_y_with_x(x, y)
#' }
.name_y_with_x <- function(x, y) {
    names(y) <- x
    return(y)
}


#' update the 'whichTree' column in row/column link data
#' @importFrom stats setNames
#' @keywords internal
#' @examples 
#' \dontrun{
#' (ld <- LinkDataFrame(nodeLab = letters[1:5],
#'                      nodeLab_alias = LETTERS[1:5],
#'                      nodeNum = 1:5,
#'                      isLeaf = TRUE,
#'                      whichTree = LETTERS[1:5],
#'                      right = 1:5))
#' newWhich <- setNames(letters[1:5], LETTERS[1:5])
#' .update_whichTree(ld, y = newWhich)
#' }
.update_whichTree <- function(x, y){
    if (is.null(x)) {
        return(x)
    }
    xx <- DataFrame(x)
    ny <- setNames(y[xx$whichTree], NULL)
    xx$whichTree <- ny
    as(xx, "LinkDataFrame")
}

#' Any element in the list is NULL
#' @keywords internal
.any_null_in_list <- function(x) {
    xl <- lapply(x, is.null)
    any(unlist(xl), na.rm = TRUE)
}


# drop tree & link data
.drop_link <- function(args, drop.colLinks, drop.rowLinks){
    # Decide whether to drop tree & link in the column dimension
    if (drop.colLinks) {
        args <- lapply(args, function(x){
            x@colTree <- NULL
            x@colLinks <- NULL
            return(x)
        })
    }
    
    # Decide whether to drop tree & link in the row dimension
    if (drop.rowLinks) {
        args <- lapply(args, function(x){
            x@rowTree <- NULL
            x@rowLinks <- NULL
            return(x)
        })
    }
    return(args)
}

#' The links & trees in the specified dim are consistent
#' @keywords internal
.is_equal_link <- function(args, dim = "row") {
    if (dim == "col") {
        link <- lapply(args, colLinks)
        tree <- lapply(args, FUN = function(x) {
            xx <- colTree(x, whichTree = NULL)})
    } else {
        link <- lapply(args, rowLinks)
        tree <- lapply(args, FUN = function(x) {
            rowTree(x, whichTree = NULL)})
    }
    
    # a list of phylo
    tree <- unlist(tree, recursive = FALSE)
    
    # all tse in args have the same tree & link in (col/row) dim
    eqL <- .all_equal_in_list(link)
    eqT <- .all_equal_in_list(tree) | is.null(tree)
    isEq <- eqT & eqL
    
    return(isEq)
}


#' rename a list automatically to avoid duplicated names
#' @keywords internal
.auto_rename_list <- function(x) {
    if (is.null(x)) { return(x)}
    names(x) <- make.names(names(x), unique = TRUE, allow_ = TRUE)
    return(x)
}



#' match a phylo to a list of phylo
#' @keywords internal
.match_phylo <- function(phy, phys) {
    ll <- lapply(phys, identical, y = phy)
    ind <- which(unlist(ll))[1]
    return(ind)
}

#' match a list of phylo (x.phys) against to a list of phylo (y.phys)
#' @keywords internal
.match_phylo_list <- function(x.phys, y.phys) {
    ll <- lapply(x.phys, .match_phylo, phys = y.phys)
    ul <- unlist(ll)
    names(ul) <- names(x.phys)
    return(ul)
}


.update_link_tree <- function(link_list, tree_list) {
    
    # new tree_list: unnest & remove duplicated trees & rename tree
    names(tree_list) <- NULL
    ntL <- unlist(tree_list, recursive = FALSE)
    oname <- names(ntL)
    ntL <- ntL[!duplicated(ntL)]
    ntL <- .auto_rename_list(x = ntL)
    
    # pair names of old & new tree_list
    ind <- lapply(tree_list, .match_phylo_list, y.phys = ntL)
    pair <- lapply(ind, FUN = function(x) {
        setNames(names(ntL)[x], names(x))
    })
    
    # update whichTree in the link data corresponding to ntL
    nlL <- mapply(.update_whichTree, link_list, pair)
    
    # new link data and list of trees
    out <- list(new_link = nlL, new_tree = ntL)
    return(out)
}

#' bind links & trees when combine TSE
#' @keywords internal
.bind_link_tree <- function(x, args, 
                            drop.rowLinks, drop.colLinks,
                            bind = "cbind") {
    
    if (bind == "rbind") { dim <- "row" } else { dim <- "col"}
    
    # Decide whether to drop tree & link
    args <- .drop_link(args = args,
                       drop.colLinks = drop.colLinks,
                       drop.rowLinks = drop.rowLinks)
    
    # old trees & links
    if (dim == "row") {
        orL <- lapply(args, rowLinks)
        otL <- lapply(args, rowTree, whichTree = NULL)
    } else {
        orL <- lapply(args, colLinks)
        otL <- lapply(args, colTree, whichTree = NULL)
    }
    
    # new trees & links (duplicated trees are removed)
    if (.all_null_in_list(otL)) {
        nT <- nL <- NULL
    } else {
        out <- .update_link_tree(link_list = orL, tree_list = otL)
        nL <- do.call(rbind, out$new_link)
        nT <- out$new_tree
    }
    
    # update slots
    if (bind == "rbind") {
        BiocGenerics:::replaceSlots(x,
                                    rowLinks = nL,
                                    rowTree = nT)
    } else {
        BiocGenerics:::replaceSlots(x,
                                    colLinks = nL,
                                    colTree = nT)
    }
    
}


#' test all TSEs have DNAStringSet in the referenceSeq slot
#' @keywords internal
.all_have_DNAStringSet <- function(args){
    refSeq <- lapply(args, FUN = function(x) {
        is(x@referenceSeq, "DNAStringSet")
    })
    all(unlist(refSeq))
}

#' test all TSEs have DNAStringSetList in the referenceSeq slot
#' @keywords internal
.all_have_DNAStringSetList <- function(args){
    refSeq <- lapply(args, FUN = function(x) {
        is(x@referenceSeq, "DNAStringSetList")
    })
    all(unlist(refSeq))
}

#' rbind referenceSeq
#' @keywords internal
#' @importFrom methods is
.rbind_refSeq <- function(args) {
    
    # all TSEs have NULL in the referenceSeq slot
    seqList <-  lapply(args, FUN = function(x) {x@referenceSeq})
    isNull <- .all_null_in_list(seqList)
    if (isNull) {return(NULL)}
    
    isDNA <- .all_have_DNAStringSet(args)
    isDNAList <- .all_have_DNAStringSetList(args)
    
    # To run rbind successfually, in the referenceSeq slot:
    #   1) all TSEs have DNAStringSet 
    #   2) all TSEs have DNAStringSetList
    #   3) all TSEs have NULL
    isV <- isNull | isDNA | isDNAList 
    
    if (!isV) {
        stop("all TSEs should have the same class in the referenceSeq slot",
             "NULL/DNAStringSet/DNAStringSetList ")
    }
    
    if (isDNA) {
        out <- do.call(c, seqList)
        return(out)
    } 
    
    if (isDNAList) {
        out <- do.call(pc, seqList)
        return(out)
    }
}

#' convert char. indicator to num. indicator
#' 
#' This differs to \code{match} with that the duplicated values in dy are not
#' ignored.
#' 
#' @param x A vector. The values to be matched.
#' @param dy A vector. The values to be matched agaist. 
#' 
#' @keywords internal
#' @author Ruizhu Huang
.match_x_dupY <- function(x, dy) {
    ul <- lapply(x, FUN = function(x) { which(dy %in% x)})
    unlist(ul)
}

#' convert char. indicator to num. indicator
#' 
#' @param ij A character or numeric indicator on rows/columns of \code{x}
#' @param x It provides row/col names for \code{ij} to be matched against.
#' @param dim "row" or "col" to specify row/col names of \code{x} to be matched
#'   against.
#' @keywords internal
#' @importFrom S4Vectors head
#' @author Ruizhu Huang
.numeric_ij <- function(ij, x, dim = "row") {
    # row/col names
    if (dim == "row") {
        char_name <- rownames(x)
    } else {
        char_name <- colnames(x)
    }
    
    if(!is.character(ij)) {return(ij)}
    
    # convert to numeric indicator
    isA <- all(ij %in% char_name)
    dff <- setdiff(ij, char_name)
    if (!isA) {
        stop(length(dff), " specified ", dim, "s can't be found.",
             call. = FALSE)
    }
    len <- sum(char_name %in% ij, na.rm = TRUE)
    ij <- match(ij, char_name)
    
    if (len > length(ij)) {
        warning("For rows/cols with the same name, only one is output")
    }
    return(ij)
}


