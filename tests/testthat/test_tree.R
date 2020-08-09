context("tree functions")

test_that("tree functions", {
    data("tinyTree")
    expect_equal(distNode(tree = tinyTree, node = c("Node_12","Node_15")),
                 0.2186453,
                 tolerance = 10^-6)
    expect_equal(distNode(tree = tinyTree, node = c("Node_13","Node_15")),
                 0.3290059,
                 tolerance = 10^-6)
    expect_equal(distNode(tree = tinyTree, node = c("Node_14","Node_15")),
                 0.6469696,
                 tolerance = 10^-6)
})
