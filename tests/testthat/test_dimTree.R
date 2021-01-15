context("rowTree & colTree")

set.seed(2)
tse_a <- makeTSE(include.colTree = FALSE)
tse_b <- makeTSE(include.colTree = FALSE)
tse_ab <- rbind(tse_a, tse_b)

tree2 <- ape::rtree(ncol(tse_a))
tree2$tip.label <- colnames(tse_a)

tree3 <- ape::rtree(nrow(tse_a))
tree3$tip.label <- rownames(tse_a)

tse_c <- makeTSE()
tse_d <- makeTSE()
tse_cd <- cbind(tse_c, tse_d)

test_that("rowTree/colTree works properly as a getter", {
    ## export one of trees
    # by the numeric index
    rt_a <- rowTree(tse_a)
    rt_b <- rowTree(tse_b)
    expect_equal(rowTree(tse_ab), rt_a) 
    expect_equal(rowTree(tse_ab, whichTree = 2), rt_b)
    
    expect_equal(colTree(tse_cd, whichTree = 1), colTree(tse_c))
    expect_equal(colTree(tse_cd, whichTree = 2), colTree(tse_d))
    ct <- colTree(tse_cd, whichTree = NULL)
    expect_equal(length(ct), 2)
    
    
    # by the name
    rnam <- rowTreeNames(tse_ab)
    expect_equal(rowTree(tse_ab, whichTree = rnam[1]), rt_a) 
    expect_equal(rowTree(tse_ab, whichTree = rnam[2]), rt_b) 
    
    
    cnam <- colTreeNames(tse_cd)
    expect_equal(colTree(tse_cd, whichTree = cnam[1]), colTree(tse_c)) 
    expect_equal(colTree(tse_cd, whichTree = cnam[2]), colTree(tse_d)) 
    
    
    ## export all trees
    rt <- rowTree(tse_ab, whichTree = NULL)
    expect_equal(length(rt), 2L)
    expect_equal(rowTree(tse_ab, whichTree = 1:2), rt)
    expect_equal(rowTree(tse_ab, whichTree = rnam), rt)
    
    ct <- colTree(tse_cd, whichTree = NULL)
    expect_equal(length(ct), 2L)
    expect_equal(colTree(tse_cd, whichTree = 1:2), ct)
    expect_equal(colTree(tse_cd, whichTree = cnam), ct)
    
 })


test_that("rowTree/colTree works properly as a setter", {
    colTree(tse_a) <- tree2
    expect_equal(colTree(tse_a), tree2)
    
    rowTree(tse_ab, whichTree = 2) <- tree3
    expect_equal(rowTree(tse_ab, whichTree = 2), tree3)
    
    
    rowTree(tse_ab, whichTree = 1) <- tree3
    expect_equal(rowTree(tse_ab, whichTree = 1), tree3)
    
    rowTree(tse_ab, whichTree = NULL) <- tree3
    expect_equal(rowTree(tse_ab, whichTree = 1), tree3)
    expect_equal(length(rowTree(tse_ab, whichTree = NULL)), 1)
})
