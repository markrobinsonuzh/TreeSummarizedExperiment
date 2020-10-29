context("MicrobiomeExperiment subset methods")

test_that("Subsetting", {
    counts <- matrix(rbinom(20, 5, .5), ncol = 4)
    rowData <- DataFrame()
    expect_error(MicrobiomeExperiment(SimpleList(counts = counts), rowData = rowData))
    rowData <- DataFrame(test = seq.int(1,5))
    refSeq <- DNAStringSetList(one = DNAStringSet(c("A","G","A","A","A")),
                               two = DNAStringSet(c("A","G","A","A","A")))
    me <- MicrobiomeExperiment(SimpleList(counts = counts), rowData = rowData,
                               referenceSeq = refSeq)

    expect_identical(me[], me)
    expect_identical(me[,], me)
    expect_identical(me[TRUE, ], me)
    expect_identical(me[, TRUE], me)
    expect_identical(me[TRUE, TRUE], me)
    actual <- me[1,]
    expect_equal(nrow(actual), 1L)
    me[1,] <- me[2,]
    expect_equal(as.character(referenceSeq(me)[[1L]][[1L]]), "G")
    me[c(1,2),] <- me[c(3,4),]
    expect_equal(as.character(referenceSeq(me)[[1L]][[1L]]), "A")
    me2 <- me
    referenceSeq(me2) <- referenceSeq(me2)[1L]
    expect_error(me2[c(1,2),] <- me[c(3,4),],
                 paste0("DNAStringSetList as 'referenceSeq' must have the ",
                        "same length to be merged"))
    referenceSeq(me2) <- NULL
    expect_error(me2[c(1,2),] <- me[c(3,4),],
                 "'x' and 'value' must have the same type of referenceSeq")
})
