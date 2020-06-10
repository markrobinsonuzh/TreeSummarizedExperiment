context("accessor_constructor")

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
tse <- TreeSummarizedExperiment(assays = list(toyTable),
                                rowData = rowInf,
                                rowNodeLab = tinyTree$tip.label,
                                colData = colInf,
                                rowTree = tinyTree)
test_that("TreeSummarizedExperiment constuctor works", {
    # TreeSummarizedExperiment
    expect_is(tse,
              "TreeSummarizedExperiment")

    # without tree
    expect_is(TreeSummarizedExperiment(assays = list(toyTable),
                                       rowData = rowInf,
                                       colData = colInf),
              "TreeSummarizedExperiment")

    expect_warning(TreeSummarizedExperiment(assays = list(toyTable),
                                          rowData = rowInf,
                                          colData = colInf,
                                          colTree = tinyTree))
})



test_that("assays could extract table successfully", {
    expect_equal(assays(tse)[[1]], toyTable)
})

test_that("assays could be written successfully", {
    assays(tse)[[2]] <- 2*toyTable
    expect_equal(assays(tse)[[2]], 2*toyTable)
})

test_that("row data could be extracted successfully", {
    expect_equal(colnames(rowData(tse)), colnames(rowInf))
    expect_setequal(colnames(rowLinks(tse)), c("nodeLab", "nodeLab_alias",
                                         "nodeNum", "isLeaf"))
    expect_equal(rowTree(tse), tinyTree)

})

test_that("row data could be written successfully", {
    tse1 <- tse
    rowData(tse1)$test <- rep(1, nrow(tse))
    expect_equal(rowData(tse1)$test, rep(1, nrow(tse)))
})


test_that("column data could be extracted successfully", {
    expect_equal(colData(tse), colInf)
})


test_that("row data could be subset successfully", {
    expect_equal(nrow(tse[1,]),1L)
    expect_equal(nrow(tse[rownames(tse)[1],]),1L)
})
test_that("column data could be subset successfully", {
    expect_equal(ncol(tse[,1]),1L)
    expect_equal(ncol(tse[,colnames(tse)[1]]),1L)
})
