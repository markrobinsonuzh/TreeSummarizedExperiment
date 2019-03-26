context("aggValue")


test_that("aggValue works correctly", {

    expect_agg_equal <- function(level, truth) {
    # assays data
    toyTable <- matrix(1:40, nrow = 10)
    colnames(toyTable) <- paste(rep(LETTERS[1:2], each = 2), rep(1:2, 2), sep = "_")
    rownames(toyTable) <- paste("entity", seq_len(10), sep = "")

    data("tinyTree")

    # row data
    rowInf <- DataFrame(var1 = sample(letters[1:2], 10, replace = TRUE),
                        var2 = sample(c(TRUE, FALSE), 10, replace = TRUE),
                        nodeLab = tinyTree$tip.label,
                        row.names = rownames(toyTable))
    # column data
    colInf <- DataFrame(gg = c(1, 2, 3, 3),
                        group = rep(LETTERS[1:2], each = 2),
                        row.names = colnames(toyTable))

    tse <- TreeSummarizedExperiment(assays = list(toyTable),
                                    rowData = rowInf,
                                    colData = colInf,
                                    rowTree = tinyTree)
    xx <- aggValue(x = tse, rowLevel = level)
    expect_setequal(assays(xx)[[1]], truth)
    }


    toyTable <- matrix(1:40, nrow = 10)

    expect_agg_equal(14, c(5, 25, 45, 65))
    expect_agg_equal(11, apply(toyTable, 2, sum))
    expect_agg_equal(c(12, 11),
                     rbind(apply(toyTable[-10, ], 2, sum),
                           apply(toyTable, 2, sum)))


})
