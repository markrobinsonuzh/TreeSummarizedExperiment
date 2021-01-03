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
    any(unlist(xl))
}

#' all elements in the list are NULL
#' @keywords internal
.all_null_in_list <- function(x) {
    xl <- lapply(x, is.null)
    all(unlist(xl))
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


#' combine trees in a list
#' @keywords internal
.cb_tree <- function(args, dim = "row", unnest = FALSE) {
    if (dim == "row") {
        out <- lapply(args, rowTree, whichTree = NULL)
    } else {
        out <- lapply(args, colTree, whichTree = NULL)
    }
    if (!unnest) {return(out)}
    
    out <- unlist(out, recursive = FALSE)
    return(out)
}

#' rename a list automatically to avoid duplicated names
#' @keywords internal
.auto_rename_list <- function(x) {
    if (is.null(x)) { return(x)}
    names(x) <- make.names(names(x), unique = TRUE, allow_ = TRUE)
    return(x)
}

#' rbind the linkDataFrame
#' @keywords internal
.cb_links <- function(args, dim = "row") {
    
    # old trees & links
    if (dim == "row") {
        orL <- lapply(args, rowLinks)
    } else {
        orL <- lapply(args, colLinks)
    }
    oldTreeList <- .cb_tree(args = args, dim = dim, unnest = FALSE)
    
    # new trees & links (duplicated trees are removed)
    if (.all_null_in_list(oldTreeList)) {return(NULL)} 
    newTrees <- unlist(oldTreeList, recursive = FALSE)
    newTrees <- newTrees[!duplicated(newTrees)]
    names(newTrees) <- make.names(names(newTrees), 
                                  unique = TRUE, allow_ = TRUE)
    
    # pair old names and new names of trees
    ind <- lapply(oldTreeList, .match_phylo_list, y.phys = newTrees)
    pair <- lapply(ind, FUN = function(x) {
        setNames(names(newTrees)[x], names(x))
     })
    
    # update the 'whichTree' column in the link data
    nrL <- mapply(.update_whichTree, orL, pair)
    
    # stack up the link data
    out <- do.call(rbind, nrL)
    return(out)
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


#' update links & trees when combine TSE
#' @keywords internal
.combine_link_tree <- function(x, args, 
                               drop.rowLinks, drop.colLinks,
                               bind = "cbind") {
    
    if (bind == "rbind") { dim <- "row" } else { dim <- "col"}
    
    # Decide whether to drop tree & link
    args <- .drop_link(args = args,
                       drop.colLinks = drop.colLinks,
                       drop.rowLinks = drop.rowLinks)
    
    
    
    # combine trees from TSEs as a list
    nT <- oT <- .cb_tree(args, dim = dim, unnest = TRUE)
    if (!is.null(oT)) {
        nT <- .auto_rename_list(x = oT)
        
        # check whether duplicated row trees exist
        if(anyDuplicated(nT)) {
            nT <- nT[!duplicated(nT)]
        }
    }
    
    # update links 
    nL <- .cb_links(args = args, dim = dim)
    
    
    
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