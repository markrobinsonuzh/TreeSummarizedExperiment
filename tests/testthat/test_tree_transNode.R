context("transNode")

data("tinyTree")

test_that("give errors when the provided arguments are not in correct form", {
    expect_error(transNode(tree = 2, node = 2))
    expect_error(transNode(tree = tinyTree, node = 20))
})

test_that("transNode could return correct results", {
    expect_transNode_equal <- function(test, truth){
        expect_equal(transNode(tree = tinyTree, node = test), truth)
    }
    expect_transNode_equal(17, "Node_17")
    expect_transNode_equal("Node_16", c(Node_16 = 16))
    expect_transNode_equal(c(2, 4, 5), c("t7", "t9", "t4"))
    expect_transNode_equal(2, "t7")
})




