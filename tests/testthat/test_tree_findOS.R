context("findOS")

data("tinyTree")

test_that("findOS return error when the node is not in the tree", {
    expect_error(findOS(node = 20, tree = tinyTree, only.leaf = TRUE))
    expect_error(findOS(node = 2, tree = 2, only.leaf = TRUE))
})

test_that("findOS could find descendant leaves correctly", {
    expect_findOS_leaf_setequal <- function(node, truth) {
        x <- findOS(node = node, tree = tinyTree,
                    only.leaf = TRUE, self.include = TRUE)[[1]]
        setequal(x, truth)
    }

    expect_findOS_leaf_setequal(15, 4:9)
    expect_findOS_leaf_setequal(13, 1:3)
    expect_findOS_leaf_setequal(18, 4:5)
    expect_findOS_leaf_setequal(11, 1:10)
})

test_that("findOS could find all descendant nodes correctly", {
    expect_findOS_node_setequal <- function(node, truth) {
        x <- findOS(node = node, tree = tinyTree,
                    only.leaf = FALSE, self.include = FALSE)[[1]]
        setequal(x, truth)
    }

    expect_findOS_node_setequal(15, c(16:19, 4:9))
    expect_findOS_node_setequal(17, c(18, 4:6))
    expect_findOS_node_setequal(11, c(1:10, 12:19))
})

test_that("findOS could include the input node in the result", {
    expect_findOS_self_setequal <- function(node, truth) {
        x <- findOS(node = node, tree = tinyTree,
                    only.leaf = TRUE, self.include = TRUE)[[1]]
        setequal(x, truth)
    }

    expect_findOS_self_setequal(5, 5)
    expect_findOS_self_setequal(17, c(17:18, 4:6))
})

