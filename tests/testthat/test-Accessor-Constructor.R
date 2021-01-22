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
colTree <- rtree(4)
colTree$tip.label <- colnames(toyTable)
colTree$node.label <- c("All", "GroupA", "GroupB")

refSeq <- DNAStringSetList(one = DNAStringSet(rep("A",nrow(toyTable))),
                           two = DNAStringSet(rep("B",nrow(toyTable))))

tse <- TreeSummarizedExperiment(assays = list(toyTable),
                                rowData = rowInf,
                                rowNodeLab = tinyTree$tip.label,
                                colData = colInf,
                                rowTree = tinyTree,
                                colTree = colTree,
                                referenceSeq = refSeq)
test_that("TreeSummarizedExperiment constuctor works", {
    # TreeSummarizedExperiment
    expect_s4_class(tse,
                    "TreeSummarizedExperiment")
    
    # without tree and refseq
    expect_s4_class(TreeSummarizedExperiment(assays = list(toyTable),
                                             rowData = rowInf,
                                             colData = colInf),
                    "TreeSummarizedExperiment")
    
    # without tree
    expect_s4_class(TreeSummarizedExperiment(assays = list(toyTable),
                                             rowData = rowInf,
                                             colData = colInf,
                                             referenceSeq = refSeq),
                    "TreeSummarizedExperiment")

    expect_warning(TreeSummarizedExperiment(assays = list(toyTable),
                                            rowData = rowInf,
                                            colData = colInf,
                                            colTree = tinyTree))
})
test_that("TreeSummarizedExperiment coercion works", {
    # SummarizedExperiment
    se <- SummarizedExperiment(assays = list(toyTable),
                               rowData = rowInf,
                               colData = colInf)
    expect_s4_class(as(se,"TreeSummarizedExperiment"),
                    "TreeSummarizedExperiment")
    # RangedSummarizedExperiment
    expect_s4_class(as(as(se,"RangedSummarizedExperiment"),
                       "TreeSummarizedExperiment"),
                    "TreeSummarizedExperiment")
    # RangedSummarizedExperiment
    expect_s4_class(as(as(se,"SingleCellExperiment"),
                       "TreeSummarizedExperiment"),
                    "TreeSummarizedExperiment")
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
                                         "nodeNum", "isLeaf", "whichTree"))
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
test_that("subsetting by node successfully", {
    expect_equal(nrow(subsetByNode(tse,"t2")),1L)
    expect_equal(subsetByNode(tse,"t2"),
                 subsetByNode(tse,1))

    expect_equal(ncol(subsetByNode(tse, colNode = "A_1")),1L)
    expect_equal(subsetByNode(tse, colNode = "A_1"),
                 subsetByNode(tse, colNode = 1))
})

test_that("other accessors/setters work", {
    # rowTreeNames & colTreeNames
    expect_equal(rowTreeNames(tse), "phylo")
    rowTreeNames(tse) <- "first_tree"
    expect_equal(rowTreeNames(tse), "first_tree")
    expect_equal(unique(rowLinks(tse)$whichTree), "first_tree")
    
    expect_equal(colTreeNames(tse), "phylo")
    colTreeNames(tse) <- "second_tree"
    expect_equal(colTreeNames(tse), "second_tree")
    expect_equal(unique(colLinks(tse)$whichTree), "second_tree")
    
    # check [
    expect_error(tse["entity11",],
                 "specified rows can't be found")
    expect_error(tse[,"C_1"],
                 "specified cols can't be found")
    
    # check show
    expect_output(show(tse),
                  "rowLinks: a LinkDataFrame")
    expect_output(show(tse),
                  "rowTree:")
    expect_output(show(tse),
                  "colLinks: a LinkDataFrame")
    expect_output(show(tse),
                  "colTree:")
    
    
    x <- TreeSummarizedExperiment(assays = list(toyTable),
                                  rowData = rowInf,
                                  colData = colInf)
    expect_output(show(x),
                  "rowLinks: NULL")
    expect_output(show(x),
                  "rowTree: NULL")
    expect_output(show(x),
                  "colLinks: NULL")
    expect_output(show(x),
                  "colTree: NULL")
    
    # rownames & colnames
    rn <- paste("entity", seq.int(10L,19L), sep = "")
    expect_true(all(rownames(rowLinks(tse)) != rn))
    rownames(tse) <- rn
    expect_true(all(rownames(rowLinks(tse)) == rn))
    cn <- paste(rep(LETTERS[1:2], each = 2), rep(3:4, 2), sep = "_")
    expect_true(all(rownames(colLinks(tse)) != cn))
    colnames(tse) <- cn
    expect_true(all(rownames(colLinks(tse)) == cn))

    tse <- TreeSummarizedExperiment(assays = list(toyTable),
                                    rowData = rowInf,
                                    colData = colInf)
    rownames(tse) <- rn
    colnames(tse) <- cn
    expect_true(all(rownames(tse) == rn))
    expect_true(all(colnames(tse) == cn))
    
    # reference seq
    expect_error(referenceSeq(tse) <- refSeq[list(1:5,1:4)])
    expect_error(referenceSeq(tse) <- refSeq[list(1:4,1:4)])
    expect_error(referenceSeq(tse) <- refSeq[[1L]][1:4])
    
    referenceSeq(tse) <- refSeq
    expect_s4_class(referenceSeq(tse),"DNAStringSetList")
    referenceSeq(tse) <- as.list(refSeq)
    expect_s4_class(referenceSeq(tse),"DNAStringSetList")
    referenceSeq(tse) <- refSeq[[1L]]
    expect_s4_class(referenceSeq(tse),"DNAStringSet")
    referenceSeq(tse) <- as.character(refSeq[[1L]])
    expect_s4_class(referenceSeq(tse),"DNAStringSet")
    
    expect_true(all(rownames(tse) == names(referenceSeq(tse))))
    rownames(tse) <- NULL
    expect_true(is.null(names(referenceSeq(tse))))
    
})
