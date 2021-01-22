context("joinode")

test_that("joinNode works properly", {
    data("tinyTree")
    expect_join_equal <- function(node, truth){
        expect_equal(joinNode(node = node, 
                              tree = tinyTree,
                              use.alias = FALSE),
                     truth)
    }
    expect_join_equal(node = c('t4','t9'), c("Node_18" = 18))
    expect_join_equal(node = c('t10','Node_18', 't8'), 
                      c("Node_17" = 17, "t10" = 7))
    expect_join_equal(node = c(5, 10), c("t4" = 5, "t3" = 10))
    expect_error(joinNode(node = c('t4','t9', 't11'), tree = tinyTree,
                           use.alias = FALSE))
    expect_warning(signalNode(node = 4:5, tree = tinyTree),
                   "'signalNode' is deprecated.")
})
