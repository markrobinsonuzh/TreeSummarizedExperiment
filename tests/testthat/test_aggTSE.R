context("aggTSE")

set.seed(1)
tse_a <- makeTSE()
tse_b <- makeTSE(include.rowTree = FALSE, include.colTree = FALSE)
tse_c <- makeTSE()
tse_ac <- rbind(tse_a, tse_c)


test_that("aggTSE works correctly", {
    # check error/warning message
    expect_warning(aggValue(x = tse_a, rowLevel = 6),
                   "'aggValue' is deprecated.")     # aggValue is deprecated
    expect_error(aggTSE(x = tse_a, rowLevel = 200)) # the node doesn't exist
    expect_error(aggTSE(x = tse_b, rowLevel = 11))  # the rowTree doesn't exist
    expect_error(aggTSE(x = tse_b, colLevel = 11))  # the rowTree doesn't exist
    expect_message(aggTSE(x = tse_a, rowLevel = 3, message = TRUE),
                   "Preparing data")
    
    # test colBlock & rowBlock
    set.seed(2)
    tse_e <- cbind(tse_a, tse_a, tse_a)
    colData(tse_e)$group <- sample(LETTERS[1:2], size = ncol(tse_e),
                                   replace = TRUE)
    
    cl <- printNode(tree = colTree(tse_e), type = "internal")$nodeNum
    tse_E <- aggTSE(x = tse_e, colLevel = cl, colBlock = "group") 
    
    expect_equal(dim(tse_E), c(10, 6))
    # For alias_7 (2 sample1 & 1 sample2 in A group)
    table(colData(tse_e)$group, colnames(tse_e))
    cd <- colData(tse_E)
    sel1 <- cd$group == "A" & colnames(tse_E) == "alias_7"
    sel2 <- cd$group == "B" & colnames(tse_E) == "alias_7"
    sel3 <- cd$group == "A" & colnames(tse_E) == "alias_5"
    sel4 <- cd$group == "B" & colnames(tse_E) == "alias_5"
    expect_equal(assays(tse_E["entity1", sel1])[[1]][1, 1], 13)
    expect_equal(assays(tse_E["entity1", sel2])[[1]][1, 1], 23)
    expect_equal(assays(tse_E["entity1", sel3])[[1]][1, 1], 117)
    expect_equal(assays(tse_E["entity1", sel4])[[1]][1, 1], 75)
    
    # aggregate on the row dim when one tree in 'rowTree'
    rlev <- printNode(tree = rowTree(tse_a), type = "all")$nodeNum
    tse_R <- aggTSE(x = tse_a, rowLevel = rlev, rowFun = sum)
    expect_equal(assays(tse_R)[[1]][14, ], 
                 setNames(c(6, 36, 66, 96), paste0("sample", 1:4)))
    
    # aggregate on the row dim when multiple trees in 'rowTree'
    rlev <- printNode(tree = rowTree(tse_ac, whichTree = 1), 
                      type = "all")$nodeNum
    tse_RR <- aggTSE(x = tse_ac, rowLevel = rlev, rowFun = sum, 
                    whichRowTree = 1)
    expect_equal(dim(tse_RR), c(19, 4))
    expect_equal(assays(tse_RR)[[1]][14, ], 
                 setNames(c(6, 36, 66, 96), paste0("sample", 1:4)))
    expect_equal(tse_RR, tse_R)
    expect_equal(colnames(rowData(tse_RR)), c("var1", "var2"))
    
    tse_RR <- aggTSE(x = tse_ac, rowLevel = rlev, rowFun = sum, 
                     whichRowTree = 1, rowDataCols = "var1")
    expect_equal(colnames(rowData(tse_RR)), "var1")
    
    # aggregate on the col dim when one tree in 'colTree'
    expect_warning(tse_ac <- cbind(tse_a, tse_c))
    clev <- printNode(tree = colTree(tse_a), type = "all")$nodeNum
    tse_C <- aggTSE(x = tse_a, colLevel = clev, colFun = sum)
    expect_equal(assays(tse_C)[[1]][, 6], 
                 setNames(seq(from = 33, by = 3, length.out = 10), 
                          paste0("entity", 1:10)))
    
    # aggregate on the col dim when multiple trees in 'colTree'
    clev <- printNode(tree = colTree(tse_ac, whichTree = 1), 
                      type = "all")$nodeNum
    tse_CC <- aggTSE(x = tse_ac, colLevel = clev, colFun = sum, 
                     whichColTree = 1)
    expect_equal(dim(tse_CC), c(10, 7))
    expect_equal(assays(tse_CC)[[1]][, 6], 
                 setNames(seq(from = 33, by = 3, length.out = 10), 
                          paste0("entity", 1:10)))
    expect_equal(tse_CC, tse_C)
    
    tse_CC <- aggTSE(x = tse_ac, colLevel = clev, colFun = sum, 
                     whichColTree = 1, colDataCols = c("ID"))
    expect_equal(colnames(colData(tse_CC)), "ID")
    
    
    clev2 <- printNode(tree = colTree(tse_ac, whichTree = 2), 
                      type = "all")$nodeNum
    tse_CC2 <- aggTSE(x = tse_ac, colLevel = clev2, colFun = sum, 
                     whichColTree = 2)
    expect_equal(dim(tse_CC2), c(10, 7))
    expect_equal(assays(tse_CC2)[[1]][, 6], 
                 setNames(seq(from = 63, by = 3, length.out = 10), 
                          paste0("entity", 1:10)))
    
    # aggregate on both dims
    rl <- printNode(tree = rowTree(tse_a, whichTree = 1), 
                    type = "internal")$nodeNum
    cl <- printNode(tree = colTree(tse_a, whichTree = 1), 
                    type = "internal")$nodeNum
    tse_B <- aggTSE(tse_a, rowLevel = rl, colLevel = cl, 
                    rowFun = sum, colFun = sum)
    count <- assays(tse_B)[[1]]
    expect_equal(count[1, 1], sum(assays(tse_a)[[1]]))
    
    # Use rowFUN = median, rowFirst = TRUE & rowFirst = FALSE
    set.seed(2)
    new_count <- matrix(rnorm(nrow(tse_a)*ncol(tse_a), mean = 10, sd = 2),
                        nrow = nrow(tse_a))
    rownames(new_count) <- rownames(tse_a)
    colnames(new_count) <- colnames(tse_a)
    assays(tse_a)[["new"]] <- new_count
    tse_RF <- aggTSE(x = tse_a, rowLevel = rl, rowFun = median, 
                     colLevel = cl, colFun = median, 
                     rowFirst = TRUE, whichAssay = "new")
    tse_CF <- aggTSE(x = tse_a, rowLevel = rl, rowFun = median, 
                     colLevel = cl, rowFirst = FALSE, whichAssay = "new")
    
    expect_false(identical(assays(tse_RF)[[1]], assays(tse_CF)[[1]]))
    
    
    
    
    # test the parallel computation
    set.seed(2)
    tse_d <- makeTSE(nrow = 1000, ncol = 1000)
    rl <- printNode(tree = rowTree(tse_d), type = "internal")$nodeNum
    cl <- printNode(tree = colTree(tse_d), type = "internal")$nodeNum
    
    tse_D <- aggTSE(x = tse_d, rowLevel = rl, colLevel = cl)
    tse_DD <- aggTSE(x = tse_d, rowLevel = rl, colLevel = cl, 
                     BPPARAM = BiocParallel::MulticoreParam(workers = 2, 
                                                            progressbar = TRUE))
    expect_equal(tse_D, tse_DD)
})



