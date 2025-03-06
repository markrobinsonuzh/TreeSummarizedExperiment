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
    xl <- lapply(x, is.null)
    xl <- lapply(xl, `!`)
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
    out <- list(new_links = nlL, new_tree = ntL)
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
        nL <- do.call(rbind, out$new_links)
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

#' replace row/col links & trees 
#' @param x A TSE with \code{ij} rows/cols to be replaced by \code{value}
#' @param value A TSE to replace some rows/cols of \code{x}.
#' @param ij A character or numeric vector to specify which rows/cols to be replaced.
#' @param dim "row" or "col" to specify the dimension is in rows or columns
#' @keywords internal
#' @author Ruizhu Huang
#' 
.replace_link_tree_1d <- function(x, value, ij, dim = "row") {
    if (missing(ij)) {
        return(NULL)
    }
    
    # multiple rows in assays might have the same name
    ij <- .numeric_ij(ij = ij, x = x, dim = dim)
    
    tseL <- list(x = x, value = value)
    if (dim == "row") {
        olL <- lapply(tseL, rowLinks)
        otL <- lapply(tseL, rowTree, whichTree = NULL)
        other <- "col"
        msg1 <- " 'rowTree()'"
        msg2 <- " 'colLinks()'"
    } else {
        olL <- lapply(tseL, colLinks)
        otL <- lapply(tseL, colTree, whichTree = NULL)
        other <- "row"
        msg1 <- " 'colTree()'"
        msg2 <- " 'rowLinks()'"
    }
    
    
    # check both w/wo tree(s) in dim
    if (!.all_null_in_list(olL) &
        !.all_nonnull_in_list(olL)) {
        stop("x' and 'value' should have the same types of", msg1,
             call. = FALSE)
    }
    # check both agrees on tree & links in the other dim
    fail_cl <- !.is_equal_link(args = tseL, dim = other)
    if (fail_cl) {
        stop("x' and 'value' differ in", msg2, call. = FALSE)
    }
    
    if (.all_nonnull_in_list(olL)) {
        # update links & trees in 'dim'
        out <- .update_link_tree(link_list = olL, 
                                 tree_list = otL)
        nlL <- out$new_links[["x"]]
        nlL[ij, ] <- out$new_links[["value"]]
        ntL <- out$new_tree[unique(nlL$whichTree)]
    } else {
        nlL <- ntL <- NULL
    }
    
    out <- list(new_links = nlL, new_tree = ntL)
    return(out)
}



.replace_link_tree_2d <- function(x, value, i, j) {
    i <- .numeric_ij(ij = i, x = x, dim = "row")
    j <- .numeric_ij(ij = j, x = x, dim = "col")
    
    tseL <- list(x = x, value = value)
    orlL <- lapply(tseL, rowLinks)
    ortL <- lapply(tseL, rowTree, whichTree = NULL)
    oclL <- lapply(tseL, colLinks)
    octL <- lapply(tseL, colTree, whichTree = NULL)
    
    
    
    # check both w/wo tree(s) in the rowLinks
    if (!.all_null_in_list(orlL) &
        !.all_nonnull_in_list(orlL)) {
        stop("x' and 'value' should have the same types of 'rowLinks()'",
             call. = FALSE)
    }
    # check both w/wo tree(s) in the colLinks
    if (!.all_null_in_list(oclL) &
        !.all_nonnull_in_list(oclL)) {
        stop("x' and 'value' should have the same types of 'colLinks()'",
             call. = FALSE)
    }
    
    # update the row link & tree
    if (.all_nonnull_in_list(orlL)) {
        # update links & trees in 'dim'
        res <- .update_link_tree(link_list = orlL, 
                                 tree_list = ortL)
        nrlL <- res$new_links[["x"]]
        nrlL[i, ] <- res$new_links[["value"]]
        nrtL <- res$new_tree[unique(nrlL$whichTree)]
    } else {
        nrlL <- nrtL <- NULL
    }
    
    
    # update the column link & tree
    if (.all_nonnull_in_list(oclL)) {
        # update links & trees in 'dim'
        res <- .update_link_tree(link_list = oclL, 
                                 tree_list = octL)
        nclL <- res$new_links[["x"]]
        nclL[j, ] <- res$new_links[["value"]]
        nctL <- res$new_tree[unique(nclL$whichTree)]
    } else {
        nclL <- nctL <- NULL
    }
    
    outR <- list(new_links = nrlL, new_tree = nrtL)
    outC <- list(new_links = nclL, new_tree = nctL)
    out <- list(outR = outR, outC = outC)
    return(out)
}

# specify which tree to be replaced In '[' replacement, i or/and j are
# specified. For example, a set of rows ('S') are mapped to a row tree ('T') .
# When 'i' is a subset of 'S', the tree ('T') can't be really removed or replace
# because there are other rows mapped to it. That is why we don't use
# .replace_link_tree_1d for the setters of rowTree/colTree

.replace_tree <- function(x, value, whichTree, 
                          nodeLab = NULL, dim = "row") {
    # Node labels of 'value'
    lab <- c(value$tip.label, value$node.label)
    empty <- c(NA, " ", "", "NA", "na")
    
    # the list of trees
    if (dim == "row") {
        tr <- rowTree(x, whichTree = NULL)
        lk <- rowLinks(x)
        nam <- rownames(x)
    } else {
        tr <- colTree(x, whichTree = NULL)
        lk <- colLinks(x)
        nam <- colnames(x)
    }
    
    # trees to be replaced
    if (!is.null(whichTree)) {
        trRep <- tr[whichTree]  
    } else {trRep <- tr}
    namRep <- names(trRep)
    
    # 'value' takes the place of the first replaced tree
    # the new list of the tree
    if (is.null(namRep[[1]])) {
        if (!is.null(tr)) {
            stop("TSE doesn't support a row/col to be mapped to multiple trees",
                 call. = FALSE)
        }
        tr <- c(tr, list("phylo" = value))
        names(tr) <- make.names(names(tr), unique = TRUE)
    } else {
        tr[[namRep[1]]] <- value
    }
    
    ntr <- tr[!names(tr) %in% namRep[-1]]
    
    # ---------------------------------------------------------------
    # update the link data
    # ---------------------------------------------------------------
    # indicate rows links to the tree to be replaced
    ii <- which(lk$whichTree %in% namRep)
    if (is.null(lk)) {
        if (dim == "row") {
            ii <- seq_len(nrow(x)) 
        } else {
            ii <- seq_len(ncol(x)) 
        }
    }
        
    if (is.null(nodeLab)) {
        olab <- nam[ii]   
    } else {olab <- nodeLab}
    
    
    # indicate rows to be dropped
    iDrop <- ii[!olab %in% lab]
    iRep <- ii[olab %in% lab]
    
    # rows has empty labels and mismatch with nodes of 'value
    mis <- olab %in% empty | !olab %in% lab
    if (sum(mis) == length(olab)) {
        stop(dim, "names of 'x' mismatch with node labels of the tree \n",
             " Try 'changeTree' with 'rowNodeLab' provided.",
             call. = FALSE)
    }
    
    if (length(iDrop)) {
        warning(length(iDrop), " ", dim, 
                "(s) are dropped due to mismatch with nodes of 'value'")
    }
    
    # update columns in the link data:
    nlk <- DataFrame(lk)
    if (!nrow(nlk)) {
        nlk <- DataFrame(
            nodeLab = olab[olab %in% lab],
            nodeNum = convertNode(tree = value, node = olab[olab %in% lab]))
        nlk$nodeLab_alias <- convertNode(tree = value, node = nlk$nodeNum, 
                                         use.alias = TRUE)
        nlk$isLeaf <- isLeaf(tree = value, node = nlk$nodeNum)
        nlk$whichTree <- names(ntr)
    } else {
        nlk$nodeLab[iRep] <- olab[olab %in% lab]
        nlk$nodeNum[iRep] <- convertNode(tree = value, node = nlk$nodeLab[iRep])
        nlk$nodeLab_alias[iRep] <- convertNode(tree = value, 
                                               node = nlk$nodeNum[iRep], 
                                               use.alias = TRUE)
        nlk$isLeaf[iRep] <- isLeaf(tree = value, 
                                   node = nlk$nodeNum[iRep])
        nlk$whichTree[iRep] <- namRep[1]
        if (length(iDrop)) {nlk <- nlk[-iDrop, ]}
    }
    
    
    # drop rows
    nlk <- as(nlk, "LinkDataFrame")
    
    out <- list(new_links = nlk, new_tree = ntr, drop = iDrop)
    
    return(out)
    
}

# This is to update old TSE objects saved using version older than 1.6.3.
# check whether the rowLinks/colLinks has the 'whichTree' column.
.lack_whichTree <- function(object, slot) {
    if (slot == "rowLinks") {lk <- rowLinks(object)}
    if (slot == "colLinks") {lk <- colLinks(object)}
    if (is.null(lk)) {return(FALSE)}
    is.null(lk$whichTree)
}

#' update dimLinks and dimTree (used in subsetByLeaf)
#' @keywords internal
#' @author Ruizhu Huang
.subset_leaf <- function(x, leaf, dim = "row", updateTree = TRUE) {
    if (dim == "row") {
        dimTree <- rowTree(x, whichTree = NULL)
        dimLink <- rowLinks(x)
    } else {
        dimTree <- colTree(x, whichTree = NULL)
        dimLink <- colLinks(x)
    }
    if (!missing(leaf)) {
        df <- lapply(seq_along(dimTree), FUN = function(ii) {
            ti <- dimTree[[ii]]
            nti <- names(dimTree)[ii]
            out <- NULL
            if (!is.numeric(leaf)) {
                lab <- intersect(leaf, c(ti$tip.label, ti$node.label))
                nd <- convertNode(tree = ti, node = lab)
            } else {
                nd <- intersect(leaf, unique(as.vector(ti$edge)))
            }
            if (length(nd)) {
                out <- data.frame(node = nd, whichTree = nti) 
            } 
            
            return(out)
        })
        df <- do.call(rbind, df)
    }
    
    ind <- which(dimLink$nodeNum %in% df$node & dimLink$whichTree %in% df$whichTree)
    if (dim == "row") { 
        x <- x[ind, ] 
        dimLink <- rowLinks(x)
        dimTree <- rowTree(x, whichTree = NULL)
    } else { 
        x <- x[, ind] 
        dimLink <- colLinks(x)
        dimTree <- colTree(x, whichTree = NULL)
    }
    
    if (!updateTree) {
        return(x)
    }
    ## update dimTree
    nam <- names(dimTree)
    new_dimTree <- lapply(nam, FUN = function(tt){
        node <- dimLink$nodeNum[dimLink$whichTree == tt]
        keep.tip(phy = dimTree[[tt]], tip = node)
    })
    names(new_dimTree) <- nam
    for (i in nam) {
        ti <- new_dimTree[[i]]
        link_i <- dimLink[dimLink[["whichTree"]] == i, "nodeLab"]
        if (dim == "row") {
            x <- changeTree(x = x, rowTree = ti, whichRowTree = i,
                            rowNodeLab = link_i)
        } else {
            x <- changeTree(x = x, colTree = ti, whichColTree = i,
                            colNodeLab = link_i)
        }
    }
    
    return(x)
}



