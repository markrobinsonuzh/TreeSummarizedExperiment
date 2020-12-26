context("TreeSummarizedExperiment subset methods")

test_that("Subsetting", {
    counts <- matrix(rbinom(20, 5, .5), ncol = 4)
    rowData <- DataFrame()
    expect_error(TreeSummarizedExperiment(SimpleList(counts = counts), rowData = rowData))
    rowData <- DataFrame(test = seq.int(1,5))
    refSeq <- DNAStringSetList(one = DNAStringSet(c("A","G","A","A","A")),
                               two = DNAStringSet(c("A","G","A","A","A")))
    tse <- TreeSummarizedExperiment(SimpleList(counts = counts), rowData = rowData,
                                    referenceSeq = refSeq)
    
    expect_identical(tse[], tse)
    expect_identical(tse[,], tse)
    expect_identical(tse[TRUE, ], tse)
    expect_identical(tse[, TRUE], tse)
    expect_identical(tse[TRUE, TRUE], tse)
    actual <- tse[1,]
    expect_equal(nrow(actual), 1L)
    tse[1,] <- tse[2,]
    expect_equal(as.character(referenceSeq(tse)[[1L]][[2L]]), "G")
    tse[c(1,2),] <- tse[c(3,4),]
    expect_equal(as.character(referenceSeq(tse)[[1L]][[1L]]), "A")
    tse2 <- tse
    referenceSeq(tse2) <- referenceSeq(tse2)[1L]
    expect_error(tse2[c(1,2),] <- tse[c(3,4),],
                 paste0("DNAStringSetList as 'referenceSeq' must have the ",
                        "same length to be merged"))
    referenceSeq(tse2) <- NULL
    expect_error(tse2[c(1,2),] <- tse[c(3,4),],
                 "'x' and 'value' must have the same type of referenceSeq")
})
