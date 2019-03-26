context("shareNode")

test_that("shareNode could find correct information", {
    data("tinyTree")
    expect_share_equal <- function(node, truth){
        expect_equal(shareNode(node = node, tree = tinyTree,
                               use.alias = FALSE),
                     truth)
    }
    expect_share_equal(node = c('t4','t9'), c("Node_18" = 18))
    expect_share_equal(node = c('t4','t9'), c("Node_18" = 18))
    expect_share_equal(node = c('t4','t9', 't1'), c("Node_16" = 16))
    expect_error(shareNode(node = c('t4','t9', 't11'), tree = tinyTree,
                           use.alias = FALSE))
})
