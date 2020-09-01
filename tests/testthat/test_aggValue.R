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
                        row.names = rownames(toyTable))
    # column data
    colInf <- DataFrame(gg = c(1, 2, 3, 3),
                        group = rep(LETTERS[1:2], each = 2),
                        row.names = colnames(toyTable))

    tse <- TreeSummarizedExperiment(assays = list(toyTable = toyTable),
                                    rowData = rowInf,
                                    colData = colInf,
                                    rowTree = tinyTree,
                                    rowNodeLab = tinyTree$tip.label)
    xx <- aggValue(x = tse, rowLevel = level)
    expect_setequal(assays(xx)[[1]], truth)
    
        expect_error(aggValue(x = tse, rowLevel = as.logical(as.numeric(as.factor(level)))),
                     "rowLevel should be a numeric or character vector.")
        expect_error(aggValue(x = tse, colLevel = as.logical(as.numeric(as.factor(level)))),
                     "colLevel should be a numeric or character vector.")
        expect_message(aggValue(x = tse, rowLevel = level, message = TRUE),
                       "Preparing data")
        expect_equal(xx, aggValue(x = tse, rowLevel = level, assay = "toyTable"))
        expect_equal(tse, aggValue(x = tse))
        
        tse2 <- TreeSummarizedExperiment(assays = list(toyTable = toyTable),
                                         rowData = rowInf,
                                         colData = colInf)
        expect_error(aggValue(x = tse2, rowLevel = level),
                     "The tree on rows doesn't exist")
        expect_error(aggValue(x = tse2, colLevel = level),
                     "The tree on columns doesn't exist")

        level_char <- transNode(tree = rowTree(tse), node = level,
                                use.alias = FALSE, message = FALSE)
        expect_equal(xx,aggValue(x = tse, rowLevel = level_char))
    }


    toyTable <- matrix(1:40, nrow = 10)

    expect_agg_equal(14, c(5, 25, 45, 65))
    expect_agg_equal(11, apply(toyTable, 2, sum))
    expect_agg_equal(c(12, 11),
                     rbind(apply(toyTable[-10, ], 2, sum),
                           apply(toyTable, 2, sum)))


})
