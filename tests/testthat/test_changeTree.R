context("changeTree")

set.seed(1)
tse_a <- makeTSE()

library(ape)
tree_1 <- tree_2 <- rtree(10)
tree_1$tip.label <- rownames(tse_a)

tree_3 <- drop.tip(tree_1, tip = rownames(tse_a)[1:3])

test_that("Test changeTree work properly", {
    # rownames of TSE can be found in node labels of the tree
    tse_a1 <- changeTree(x = tse_a, rowTree = tree_1)
    expect_s4_class(tse_a1,
                    "TreeSummarizedExperiment")
    expect_equal(nrow(tse_a1), nrow(tse_a))
    
    # rownames of TSE can't be found in node labels of the tree
    expect_error(changeTree(x = tse_a, rowTree = tree_2))
    tse_a2 <- changeTree(x = tse_a, rowTree = tree_2,
                         rowNodeLab = tree_2$tip.label)
    expect_equal(nrow(tse_a2), nrow(tse_a))
    
    # duplicated rownames; rownames agree with node labels of the tree
    tse_aa <- rbind(tse_a, tse_a)
    tse_aa1 <- changeTree(x = tse_aa, rowTree = tree_1)
    expect_equal(nrow(tse_aa1), nrow(tse_aa))
    
    # duplicated rownames; rownames can't be found in node labels of the tree
    expect_warning(tse_aa3 <- changeTree(x = tse_aa, rowTree = tree_3))
    expect_equal(nrow(tse_aa3), nrow(tse_aa) - 6)
    
})
