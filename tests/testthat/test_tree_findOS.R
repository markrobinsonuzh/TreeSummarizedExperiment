context("findOS")

data("tinyTree")

test_that("findOS return error when the node is not in the tree", {
    expect_error(findOS(node = 20, tree = tinyTree, only.leaf = TRUE))
    expect_error(findOS(node = 2, tree = 2, only.leaf = TRUE))
})

test_that("findOS could find descendant leaves correctly", {
    expect_findOS_leaf_expect_equal <- function(node, truth) {
        x <- findOS(node = node, tree = tinyTree,
                    only.leaf = TRUE, self.include = TRUE)[[1]]
        expect_equal(x, truth)
    }

    expect_findOS_leaf_expect_equal(15, 4:9)
    expect_findOS_leaf_expect_equal(13, 1:3)
    expect_findOS_leaf_expect_equal(18, 4:5)
    expect_findOS_leaf_expect_equal(11, 1:10)
})

test_that("findOS could find all descendant nodes correctly", {
    expect_findOS_leaf_expect_equal <- function(node, truth) {
        x <- findOS(node = node, tree = tinyTree,
                    only.leaf = FALSE, self.include = FALSE)[[1]]
        expect_equal(x, truth)
    }

    expect_findOS_leaf_expect_equal(15, c(4:9,16:19))
    expect_findOS_leaf_expect_equal(17, c(4:6,18))
    expect_findOS_leaf_expect_equal(11, c(1:10, 12:19))
})

test_that("findOS could include the input node in the result", {
    expect_findOS_leaf_expect_equal <- function(node, truth) {
        x <- findOS(node = node, tree = tinyTree,
                    only.leaf = TRUE, self.include = TRUE)[[1]]
        expect_equal(x, truth)
    }

    expect_findOS_leaf_expect_equal(5, 5)
    expect_findOS_leaf_expect_equal(17, c(4:6))
})

