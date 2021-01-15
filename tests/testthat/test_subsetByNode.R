context("subsetByNode")

set.seed(1)
tse_a <- makeTSE(include.rowTree = FALSE)
tse_b <- makeTSE(include.rowTree = FALSE)
tse_ab <- cbind(tse_a, tse_b)

   
test_that("subsetByNode works properly", {
    
    # the row dim.
    expect_warning(subsetByNode(x = tse_ab, colNode = 2:3,
                                whichRowTree = "phylo"),
                   "The row tree is not available")
    
    # get data linked to the tree 'phylo'
    sse_ab <- subsetByNode(x = tse_ab, colNode = 2:3,
                           whichColTree = "phylo")
    expect_equal(unique(colLinks(sse_ab)$whichTree), "phylo")
    
    # get the data linked nodes (2:3) of all trees
    sse_ab <- subsetByNode(x = tse_ab, colNode = 2:3)
    expect_equal(unique(colLinks(sse_ab)$whichTree), c("phylo", "phylo.1"))
    
    # by row & col trees
    tree2 <- ape::rtree(10)
    tree2$tip.label <- rownames(tse_ab)
    rowTree(tse_ab, whichTree = NULL) <- tree2
    
    sse <- subsetByNode(x = tse_ab, whichRowTree = 1, whichColTree = 1)
    expect_equal(dim(sse), c(10, 4))
    
    sse <- subsetByNode(x = tse_ab, rowNode = 1:5, colNode = 1:2,
                        whichRowTree = 1, whichColTree = 1)
    expect_equal(dim(sse), c(5, 2))
    
})

