context("TreeSummarizedExperiment subset methods")

test_that("Subsetting", {
    counts <- matrix(rbinom(20, 5, .5), ncol = 4)
    rowData <- DataFrame()
    expect_error(TreeSummarizedExperiment(SimpleList(counts = counts), rowData = rowData))
    rowData <- DataFrame(test = seq.int(1,5))
    refSeq <- DNAStringSetList(one = DNAStringSet(c("A","G","A","A","A")),
                               two = DNAStringSet(c("A","G","A","A","A")))
    tse <- TreeSummarizedExperiment(SimpleList(counts = counts), 
                                    rowData = rowData,
                                    referenceSeq = refSeq)
    
    expect_identical(tse[], tse)
    expect_identical(tse[,], tse)
    expect_identical(tse[TRUE, ], tse)
    expect_identical(tse[, TRUE], tse)
    expect_identical(tse[TRUE, TRUE], tse)
    expect_identical(nrow(tse[-1, ]), 4L)
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


set.seed(1)
tse_a <- makeTSE(include.colTree = FALSE)
tse_b <- makeTSE(include.colTree = FALSE)
tse_c <- makeTSE(include.colTree = FALSE)
test_that("Subset & replace works successfully in the row dim.", {
   # 'value' and 'x' has the same row tree
    tse1 <- tse_a
    tse1[1, ] <- tse1[2, ]
    expect_equal(nrow(tse1), nrow(tse_a))
    expect_equal(length(rowTree(tse1, whichTree = NULL)), 1)
    
    # 'value' and 'x' has different row trees
    nr <- nrow(tse_a)
    tse_a[1:2, ] <- tse_b[1:2, ]
    rt <- rowTree(tse_a, whichTree = NULL)
    rl <- rowLinks(tse_a)
    expect_equal(nrow(tse_a), nr)
    expect_equal(length(rt), 2)
    expect_equal(length(unique(rl$whichTree)), 2)
    
    # 'x' has multiple row trees
    tse_ab <- rbind(tse_a, tse_b)
    nr <- nrow(tse_ab)
    tse_ab[1:3, ] <- tse_c[1:3, ] 
    rt <- rowTree(tse_ab, whichTree = NULL)
    rl <- rowLinks(tse_ab)
    expect_equal(nrow(tse_ab), nr)
    expect_equal(length(rt), 3)
    expect_equal(length(unique(rl$whichTree)), 3)
    
})


set.seed(1)
tse_a <- makeTSE(include.colTree = FALSE)
tse_b <- makeTSE(include.colTree = FALSE)
tse_c <- makeTSE(include.colTree = FALSE)
test_that("Subset & replace works successfully in the row dim.", {
    # 'value' and 'x' has the same row tree
    tse1 <- tse_a
    tse1[1, ] <- tse1[2, ]
    expect_equal(nrow(tse1), nrow(tse_a))
    expect_equal(length(rowTree(tse1, whichTree = NULL)), 1)
    
    # 'value' and 'x' has different row trees
    nr <- nrow(tse_a)
    tse_a[1:2, ] <- tse_b[1:2, ]
    rt <- rowTree(tse_a, whichTree = NULL)
    rl <- rowLinks(tse_a)
    expect_equal(nrow(tse_a), nr)
    expect_equal(length(rt), 2)
    expect_equal(length(unique(rl$whichTree)), 2)
    
    # 'x' has multiple row trees
    tse_ab <- rbind(tse_a, tse_b)
    nr <- nrow(tse_ab)
    tse_ab[1:3, ] <- tse_c[1:3, ] 
    rt <- rowTree(tse_ab, whichTree = NULL)
    rl <- rowLinks(tse_ab)
    expect_equal(nrow(tse_ab), nr)
    expect_equal(length(rt), 3)
    expect_equal(length(unique(rl$whichTree)), 3)
    expect_equal(unique(rl$whichTree), rowTreeNames(tse_ab))
})



set.seed(1)
tse_a <- makeTSE(include.rowTree = FALSE)
tse_b <- makeTSE(include.rowTree = FALSE)
tse_c <- makeTSE(include.rowTree = FALSE)
tse_F <- makeTSE(nrow = 20, include.rowTree = FALSE)
test_that("Subset & replace works successfully in the col dim.", {
    # 'value' and 'x' has the same row tree
    tse1 <- tse_a
    tse1[, 1] <- tse1[, 2]
    expect_equal(ncol(tse1), ncol(tse_a))
    expect_equal(length(colTree(tse1, whichTree = NULL)), 1)
    
    # 'value' and 'x' has different col trees
    nr <- ncol(tse_a)
    tse_a[, 1:2] <- tse_b[, 1:2]
    rt <- colTree(tse_a, whichTree = NULL)
    rl <- colLinks(tse_a)
    expect_equal(ncol(tse_a), nr)
    expect_equal(length(rt), 2)
    expect_equal(length(unique(rl$whichTree)), 2)
    
    # 'x' has multiple col trees
    tse_ab <- cbind(tse_a, tse_b)
    nr <- ncol(tse_ab)
    tse_ab[, 1:3] <- tse_c[, 1:3] 
    rt <- colTree(tse_ab, whichTree = NULL)
    rl <- colLinks(tse_ab)
    expect_equal(ncol(tse_ab), nr)
    expect_equal(length(rt), 3)
    expect_equal(length(unique(rl$whichTree)), 3)
    expect_equal(unique(rl$whichTree), colTreeNames(tse_ab))
    
    
    expect_error(tse_a[, 2] <- tse_F[, 3], "number of items to replace")
})


set.seed(1)
tse_a <- makeTSE()
tse_b <- makeTSE()

test_that("Subset & replace works successfully in both row & col dim.", {
    # 'value' and 'x' has the same row & col trees 
    tse1 <- tse_a
    tse1[1, 2] <- tse_a[3, 2]
    expect_equal(dim(tse1), dim(tse_a))
    expect_equal(length(colTree(tse1, whichTree = NULL)), 1)
    expect_equal(length(rowTree(tse1, whichTree = NULL)), 1)
    expect_equal(rownames(tse1), c(rownames(tse_a)[3], rownames(tse1)[-1]))
    cname <- c(colnames(tse1)[1], colnames(tse_a)[2], colnames(tse1)[-(1:2)])
    expect_equal(colnames(tse1), cname)
    
    # 'value' and 'x' has different row trees & col trees 
    tse1 <- tse_a
    tse1[1, 2] <- tse_b[3, 2]
    expect_equal(dim(tse1), dim(tse_a))
    expect_equal(length(colTree(tse1, whichTree = NULL)), 2)
    expect_equal(length(rowTree(tse1, whichTree = NULL)), 2)
    expect_equal(rownames(tse1), c(rownames(tse_b)[3], rownames(tse1)[-1]))
    cname <- c(colnames(tse1)[1], colnames(tse_b)[2], colnames(tse1)[-(1:2)])
    expect_equal(colnames(tse1), cname)
    
    namRx <- rowLinks(tse1)$whichTree[1]
    namCx <- colLinks(tse1)$whichTree[2]
    namRv <- rowLinks(tse_b)$whichTree[3]
    namCv <- colLinks(tse_b)$whichTree[2]
    expect_identical(rowTree(tse1, whichTree = namRx),
                     rowTree(tse_b, whichTree = namRv))
    expect_identical(colTree(tse1, whichTree = namCx),
                     colTree(tse_b, whichTree = namCv))
    
    
    
})
