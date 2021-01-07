context("rbind-cbind")

# rbind works :
# a) TSE without rowTree and without colTree
# b) TSE with rowTree but without colTree
# c) TSE without rowTree but with colTree
# d) TSE with rowTree & colTree
# e) TSE with rowTree but without colTree (but has a different rowTree as b)
# f) TSE without rowTree but with colTree (but has a different colTree as c)

set.seed(1)
# a)
(tse_a <- makeTSE(include.rowTree = FALSE, include.colTree = FALSE))

# b)
(tse_b <- makeTSE(include.colTree = FALSE))
(tse_e <- makeTSE(include.colTree = FALSE))

# c)
(tse_c <- makeTSE(include.rowTree = FALSE))
(tse_f <- makeTSE(include.rowTree = FALSE))

# d)
(tse_d <- makeTSE())


test_that("cbind TreeSummarizedExperiment successfully", {
    # TreeSummarizedExperiment
    expect_s4_class(cbind(tse_c, tse_f),
                    "TreeSummarizedExperiment")
    
    # drop colTree when some TSEs have NULL and others have non-NULL colTree
    expect_warning(cbind(tse_a, tse_c), 
                   "colTree should be all NULL or non-NULL")
    
    # drop rowTree & rowLinks when rowTree are different in TSEs
    expect_warning(cbind(tse_b, tse_e), 
                   "rowTree & rowLinks differ in the provided TSEs")
    expect_warning(cbind(tse_b, tse_c), 
                   "rowTree & rowLinks differ in the provided TSEs")
    expect_warning(cbind(tse_b, tse_c), 
                   "colTree should be all NULL or non-NULL in TSEs")
    
    # cbind TSEs with different colTree end up with two phylo objects in colTree
    # of the new TSE
    tse_cf <- cbind(tse_c, tse_f)
    expect_equal(length(colTree(tse_cf, whichTree = NULL)), 2)
    
    # rbind TSEs with multiple phylo in the rowTree
    tse_cfc <- cbind(tse_cf, tse_c)
    expect_equal(length(colTree(tse_cfc, whichTree = NULL)), 2)
    expect_equal(nrow(tse_cfc), 10)
    expect_equal(ncol(tse_cfc), 12)
    
    # cbind multiple TSEs
    tse_CFC <- cbind(tse_c, tse_f, tse_c)
    expect_equal(tse_CFC, tse_cfc)
    
    # duplicated colTree are removed
    tse_cc <- cbind(tse_c, tse_c)
    expect_equal(length(colTree(tse_cc, whichTree = NULL)), 1)
})




test_that("rbind TreeSummarizedExperiment successfully", {
    # TreeSummarizedExperiment
    expect_s4_class(rbind(tse_b, tse_e),
                    "TreeSummarizedExperiment")
    
    # drop rowTree when some of TSEs have NULL and others have non-NULL rowTree
    expect_warning(rbind(tse_a, tse_b), 
                   "rowTree should be all NULL or non-NULL")
    
    # drop colTree & colLinks when colTree are different in TSEs
    expect_warning(rbind(tse_a, tse_c), 
                   "colTree & colLinks differ in the provided TSEs")
    expect_warning(rbind(tse_a, tse_d), 
                   "colTree & colLinks differ in the provided TSEs")
    expect_warning(rbind(tse_a, tse_d), 
                   "rowTree should be all NULL or non-NULL in TSEs")
    
    # rbind TSEs with different rowTree end up with two phylo objects in rowTree
    # of the new TSE
    tse_be <- rbind(tse_b, tse_e)
    expect_equal(length(rowTree(tse_be, whichTree = NULL)), 2)
    
    # rbind TSEs with multiple phylo in the rowTree
    tse_bbe <- rbind(tse_be, tse_b)
    expect_equal(length(rowTree(tse_bbe, whichTree = NULL)), 2)
    expect_equal(nrow(tse_bbe), 30)
    expect_equal(ncol(tse_bbe), 4)
    expect_equal(rowLinks(tse_bbe)$whichTree[1], "phylo")
    expect_equal(rowLinks(tse_bbe)$whichTree[11], "phylo.1")
    expect_equal(rowLinks(tse_bbe)$whichTree[21], "phylo")
    
    # rbind multiple TSEs
    tse_BBE <- rbind(tse_b, tse_e, tse_b)
    expect_equal(tse_BBE, tse_bbe)
    
    # duplicated rowTree are removed
    tse_bb <- rbind(tse_b, tse_b)
    expect_equal(length(rowTree(tse_bb, whichTree = NULL)), 1)
})

test_that("rbind TreeSummarizedExperiment with referenceSeq successfully", {
    library(BiocGenerics)
    refSeq1 <- DNAStringSetList(one = DNAStringSet(rep("A", 10)),
                                two = DNAStringSet(rep("B", 10)))
    refSeq2 <- DNAStringSetList(one = DNAStringSet(rep("C", 10)),
                                two = DNAStringSet(rep("D", 10)))
    refSeq3 <- DNAStringSet(rep("C", 10))
    
    referenceSeq(tse_b) <- refSeq1
    referenceSeq(tse_e) <- refSeq2
    
    # rbind TSEs with DNAStringSetList in the  referenceSeq slot
    tse_be <- rbind(tse_b, tse_e)
    expect_equal(lengths(referenceSeq(tse_be)), c(one = 20, two = 20))
    
    # Error: one TSE with DNAStringSetList & the other with DNAStringSet
    referenceSeq(tse_e) <- refSeq3
    expect_error(rbind(tse_b, tse_e),
                 "all TSEs should have the same class in the referenceSeq")
    
    # Error: one TSE with DNAStringSetList & the other with NULL
    referenceSeq(tse_e) <- NULL
    expect_error(rbind(tse_b, tse_e),
                 "all TSEs should have the same class in the referenceSeq")
})