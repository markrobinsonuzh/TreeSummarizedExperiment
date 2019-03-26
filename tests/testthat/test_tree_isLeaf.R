context("isLeaf")

test_that("shareNode could find correct information", {
    data("tinyTree")
    expect_isLeaf_equal <- function(test, truth){
        expect_equal(isLeaf(tree = tinyTree, node = test), truth)
    }
    expect_isLeaf_equal(5, TRUE)
    expect_isLeaf_equal(4, TRUE)
    expect_isLeaf_equal(18, FALSE)
    expect_error(isLeaf(tree = tinyTree, input = 20))
})
