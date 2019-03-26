context("findSibling")

test_that("findSibling could find correct nodes", {
     data(tinyTree)

    expect_sibling_setequal <- function(node, truth) {
        expect_setequal(findSibling(tree = tinyTree, node = node,
                                 use.alias = FALSE), truth)
    }

    expect_sibling_setequal(17, 19)
    expect_sibling_setequal(16, 9)
    expect_sibling_setequal(c(17, 16), c(19, 9))
    expect_error(findSibling(tree = tinyTree, node = 51))

})
