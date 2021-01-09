context("setter: rowTree & colTree")

set.seed(1)
n <- 10
tse_a <- makeTSE(nrow = n, include.colTree = FALSE)
tse_b <- makeTSE(nrow = n, include.colTree = FALSE)
tse_ab <- rbind(tse_a, tse_b)


new_tree1 <- ape::rtree(n)
new_tree1$tip.label <- rownames(tse_a)

new_tree2 <- ape::rtree(n/2)
new_tree2$tip.label <- rownames(tse_a)[seq_len(n/2)]

test_that("Repace the row tree successfully", {
    # Only one tree in the slot rowTree
    tse_x <- tse_a
    rowTree(x = tse_x, whichTree = 1) <- new_tree1
    expect_equal(rowTree(tse_x, whichTree = 1), new_tree1)
    expect_false(identical(rowTree(tse_x, whichTree = 1),
                           rowTree(tse_a, whichTree = 1)))
    expect_warning(rowTree(x = tse_x, whichTree = 1) <- new_tree2)
    
    
    # Two trees in the slot rowTree: replace the first tree.
    # Entity10 - 1(node), Entity9 - 2 (node), ... , Entity1 - 10(node)
    new_tree1$tip.label <- paste0("entity", 10:1)
    rtn <- rowTreeNames(tse_ab)
    rowTree(x = tse_ab, whichTree = rtn[1]) <- new_tree1
    rt <- rowTree(tse_ab, whichTree = NULL)
    rl <- rowLinks(tse_ab)
    
    
    expect_equal(length(rt), 2)
    expect_equal(rt[[1]], new_tree1)
    expect_equal(rl$nodeLab[rl$whichTree == rtn[1] & rl$nodeNum == 10],
                 "entity1")
    
    # Two trees in the slot rowTree: replace all tree with 'value'.
    rowTree(x = tse_ab, whichTree = NULL) <- new_tree1
    rt <- rowTree(tse_ab, whichTree = NULL)
    rl <- rowLinks(tse_ab)
    sub_rl <- rl[rl$nodeNum == 10, ]
    
    expect_equal(length(rt), 1)
    expect_equal(rt[[1]], new_tree1)
    expect_equal(nrow(sub_rl), 2)
    
    })


set.seed(1)
n <- 10
tse_a <- makeTSE(ncol = n, include.rowTree = FALSE)
tse_b <- makeTSE(ncol = n, include.rowTree = FALSE)
tse_ab <- cbind(tse_a, tse_b)


new_tree1 <- ape::rtree(n)
new_tree1$tip.label <- colnames(tse_a)

new_tree2 <- ape::rtree(n/2)
new_tree2$tip.label <- colnames(tse_a)[seq_len(n/2)]

test_that("Repace the column tree successfully", {
    # Only one tree in the slot rowTree
    tse_x <- tse_a
    colTree(x = tse_x, whichTree = 1) <- new_tree1
    expect_equal(colTree(tse_x, whichTree = 1), new_tree1)
    expect_false(identical(colTree(tse_x, whichTree = 1),
                           colTree(tse_a, whichTree = 1)))
    
    expect_warning(colTree(x = tse_x, whichTree = 1) <- new_tree2)
    
    
    # Two trees in the slot colTree: replace the first tree.
    # Sample10 - 1(node), Sample9 - 2 (node), ... , Sample1 - 10(node)
    new_tree1$tip.label <- paste0("sample", 10:1)
    ctn <- colTreeNames(tse_ab)
    colTree(x = tse_ab, whichTree = ctn[1]) <- new_tree1
    ct <- colTree(tse_ab, whichTree = NULL)
    cl <- colLinks(tse_ab)
    
    
    expect_equal(length(ct), 2)
    expect_equal(ct[[1]], new_tree1)
    expect_equal(cl$nodeLab[cl$whichTree == ctn[1] & cl$nodeNum == 10],
                 "sample1")
    
    # Two trees in the slot rowTree: replace all tree with 'value'.
    colTree(x = tse_ab, whichTree = NULL) <- new_tree1
    ct <- colTree(tse_ab, whichTree = NULL)
    cl <- colLinks(tse_ab)
    sub_cl <- cl[cl$nodeNum == 10, ]
    
    expect_equal(length(ct), 1)
    expect_equal(ct[[1]], new_tree1)
    expect_equal(nrow(sub_cl), 2)
    
})

