context("matTree")

test_that("findSibling could find correct nodes", {
    data(tinyTree)

    mt <- matTree(tree = tinyTree)

    expect_matTree_setequal <- function(ith, truth) {
        expect_setequal(matTree(tree = tinyTree)[ith, ], truth)
    }

    expect_matTree_setequal(1, c(1, 11:13, rep(NA, 3)))
    expect_matTree_setequal(2, c(2, 11:14, rep(NA, 2)))
    expect_matTree_setequal(3, c(3, 11:14, rep(NA, 2)))
    expect_matTree_setequal(4, c(4, 11, 12, 15:18))
    expect_matTree_setequal(10, c(10:11, rep(NA, 5)))

})
