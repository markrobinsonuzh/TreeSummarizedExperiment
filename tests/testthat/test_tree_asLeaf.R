context("asLeaf")

test_that("asLeaf", {
    data("tinyTree")
    
    expect_error(asLeaf(tinyTree),
                 'argument "node" is missing')
    expect_error(asLeaf(tree = tinyTree, node = "t2"),
                 "Leaf node\\(s\\) found")
    expect_error(asLeaf(tree = tinyTree, node = FALSE),
                 "The input node should be a character or numeric vector")
    expect_error(asLeaf(tree = tinyTree, node = "Node_11"),
                 "Selected root node of tree")
    
    expect_asLeaf_equal <- function(test, truth, truth2, truth3){
        actual <- asLeaf(tree = tinyTree, node = test)
        expect_equal(nrow(actual$edge), truth)
        expect_equal(countLeaf(actual), truth2)
        expect_equal(countNode(actual), truth3)
    }
    expect_asLeaf_equal("Node_12", 2L, 2L, 3L)
    expect_asLeaf_equal("Node_13", 14L, 8L, 15L)
    expect_asLeaf_equal("Node_14", 16L, 9L, 17L)
    expect_asLeaf_equal("Node_15", 8L, 5L, 9L)
    expect_asLeaf_equal("Node_16", 10L, 6L, 11L)
    expect_asLeaf_equal("Node_17", 14L, 8L, 15L)
    expect_asLeaf_equal("Node_18", 16L, 9L, 17L)
    expect_asLeaf_equal("Node_19", 16L, 9L, 17L)
    expect_asLeaf_equal(c("Node_12","Node_15"), 2L, 2L, 3L)
    expect_asLeaf_equal(c("Node_13","Node_15"), 4L, 3L, 5L)
    expect_asLeaf_equal(c("Node_14","Node_15"), 6L, 4L, 7L)
    expect_equal(asLeaf(tree = tinyTree, node = c("Node_12","Node_15")),
                 asLeaf(tree = tinyTree, node = c("Node_15","Node_12")))
    expect_error(asLeaf(tree = tinyTree, node = "Node_20"),
                 "1 nodes mismatch with the tree")
})

