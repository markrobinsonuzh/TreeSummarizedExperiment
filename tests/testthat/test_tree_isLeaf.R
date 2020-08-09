context("isLeaf")

test_that("shareNode could find correct information", {
    data("tinyTree")
    expect_error(isLeaf(tree = tinyTree),
                 'argument "node" is missing')
    
    expect_isLeaf_equal <- function(test, truth){
        expect_equal(isLeaf(tree = tinyTree, node = test), truth)
    }
    expect_isLeaf_equal(5, TRUE)
    expect_isLeaf_equal(4, TRUE)
    expect_isLeaf_equal(18, FALSE)
    expect_isLeaf_equal(c(5,4,18), c(TRUE,TRUE,FALSE))
    expect_isLeaf_equal(c(5,18,4), c(TRUE,FALSE,TRUE))
    expect_error(isLeaf(tree = tinyTree, node = 20),
                 "Node20 can't be matched to any node of the tree")
    expect_error(isLeaf(tree = tinyTree, input = 20),
                 "unused argument \\(input = 20\\)")
})
