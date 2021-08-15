context("subsetByLeaf")
library(ape)
set.seed(1)
z <- makeTSE(nrow = 5, ncol = 4, include.rowTree = TRUE, include.colTree = FALSE)
y <- makeTSE(nrow = 4, ncol = 4, include.rowTree = TRUE, include.colTree = FALSE)
tr <- ape::rtree(4)
zy <- rbind(z, y)

tse <- changeTree(x = zy, rowTree = tr, whichRowTree = 2, rowNodeLab = tr$tip.label)
   
test_that("subsetByLeaf works properly in the row dimension", {
    
    # warning if the dim tree is not available.
    expect_warning(subsetByLeaf(x = tse, colLeaf = 2:3,
                                whichRowTree = "phylo"),
                   "The column tree is not available")
    
    ## 1) rowLeaf exist only in one of trees
    rf <- c("t4", "t2")
    sx <- subsetByLeaf(x = tse, rowLeaf = rf, updateTree = TRUE)
    expect_equal(rowLinks(sx)$nodeLab, rf)
    expect_equal(unname(rowLinks(sx)$nodeNum), c(1, 2))
    
    sx <- subsetByLeaf(x = tse, rowLeaf = rf, updateTree = FALSE)
    expect_equal(unname(rowLinks(sx)$nodeNum), c(1, 4))
    
    ## 2) rowLeaf exist in all trees
    rf <- 1:3
    sx <- subsetByLeaf(x = tse, rowLeaf = rf)
    expect_equal(length(rowTree(x = sx, whichTree = NULL)), 2)
    expect_equal(rowLinks(sx)$nodeNum, 
                 setNames(rep(1:3, 2), paste0("entity", rep(1:3, 2))))
    
    
    ## 3) rowLeaf exist in all trees, but subset and update 
    ##    only the specified
    rf <- 1:3
    sx <- subsetByLeaf(x = tse, rowLeaf = rf, whichRowTree = "phylo")
    expect_equal(length(rowTree(x = sx, whichTree = NULL)), 1)
    expect_equal(rowLinks(sx)$nodeNum, 
                 setNames(1:3, paste0("entity", 1:3)))
    
    
    rf <- 3:4
    sx <- subsetByLeaf(x = tse, rowLeaf = rf, whichRowTree = "phylo.1")
    expect_equal(length(rowTree(x = sx, whichTree = NULL)), 1)
    expect_equal(rowLinks(sx)$nodeNum, 
                 setNames(1:2, paste0("entity", 3:4)))
    
    rf <- 3:4
    sx <- subsetByLeaf(x = tse, rowLeaf = rf, whichRowTree = "phylo.1",
                       updateTree = FALSE)
    expect_equal(length(rowTree(x = sx, whichTree = NULL)), 1)
    expect_equal(rowLinks(sx)$nodeNum, 
                 setNames(3:4, paste0("entity", 3:4)))
})

z <- makeTSE(nrow = 5, ncol = 5, include.rowTree = FALSE, include.colTree = TRUE)
y <- makeTSE(nrow = 5, ncol = 4, include.rowTree = FALSE, include.colTree = TRUE)
tr <- ape::rtree(4)
zy <- cbind(z, y)

tse <- changeTree(x = zy, colTree = tr, whichColTree = 2, colNodeLab = tr$tip.label)

test_that("subsetByLeaf works properly in the column dimension", {
    
    # warning if the dim tree is not available.
    expect_warning(subsetByLeaf(x = tse, rowLeaf = 2:3,
                                whichRowTree = "phylo"),
                   "The row tree is not available")
    
    ## 1) rowLeaf exist only in one of trees
    rf <- c("t2", "t4")
    sx <- subsetByLeaf(x = tse, colLeaf = rf, updateTree = TRUE)
    expect_equal(colLinks(sx)$nodeLab, rf)
    expect_equal(unname(colLinks(sx)$nodeNum), c(1, 2))
    
    sx <- subsetByLeaf(x = tse, colLeaf = rf, updateTree = FALSE)
    expect_equal(unique(colLinks(sx)$whichTree), "phylo.1")
    
    ## 2) rowLeaf exist in all trees
    rf <- 1:3
    sx <- subsetByLeaf(x = tse, colLeaf = rf)
    expect_equal(length(colTree(x = sx, whichTree = NULL)), 2)
    expect_equal(unique(colLinks(sx)$whichTree), c("phylo", "phylo.1"))
    
    
    ## 3) rowLeaf exist in all trees, but subset and update 
    ##    only the specified
    rf <- 1:3
    sx <- subsetByLeaf(x = tse, colLeaf = rf, whichColTree = "phylo")
    expect_equal(length(colTree(x = sx, whichTree = NULL)), 1)
    expect_equal(unique(colLinks(sx)$whichTree), "phylo")
})


tse <- makeTSE(nrow = 10, ncol = 5)
test_that("subsetByLeaf works properly in both dimensions", {
    
    sse <- subsetByLeaf(x = tse, rowLeaf = 1:3, colLeaf = 3:5)
    expect_equal(dim(sse), c(3, 3))
})
