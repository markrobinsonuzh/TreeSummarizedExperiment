context("tree functions")

test_that("tree functions", {
    data("tinyTree")
    # distNode
    expect_equal(distNode(tree = tinyTree, node = c("Node_12","Node_15")),
                 0.2186453,
                 tolerance = 10^-6)
    expect_equal(distNode(tree = tinyTree, node = c("Node_13","Node_15")),
                 0.3290059,
                 tolerance = 10^-6)
    expect_equal(distNode(tree = tinyTree, node = c("Node_14","Node_15")),
                 0.6469696,
                 tolerance = 10^-6)
    
    # findChild
    expect_equivalent(findChild(tree = tinyTree, node = c("Node_12","Node_15")),
                      c(15,16))
    expect_equivalent(findChild(tree = tinyTree, node = c("Node_13","Node_15")),
                      c(1,16))
    expect_equivalent(findChild(tree = tinyTree, node = c("Node_14","Node_15")),
                      c(2,16))
    
    # printNode
    expect_error(printNode(tree = "a"),
                 "tree should be a phylo object")
    expect_error(printNode(tree = tinyTree, type = "node"),
                 "'arg' should be one of ")
    check_printNode <- function(type){
        actual <- printNode(tree = tinyTree, type = "leaf")
        expect_s3_class(actual,"data.frame")
        expect_named(actual,c("nodeLab", "nodeLab_alias", "nodeNum", "isLeaf"))
    }
    check_printNode("leaf")
    check_printNode("internal")
    check_printNode("all")
    
    # showNode
    actual <- showNode(tree = tinyTree)
    expect_equal(unname(actual), seq_len(19L))
    expect_named(actual, c(tinyTree$tip.label,tinyTree$node.label))
    actual <- showNode(tree = tinyTree, only.leaf = TRUE)
    expect_equal(unname(actual), seq_len(10L))
    expect_named(actual, c(tinyTree$tip.label))
    actual <- showNode(tree = tinyTree, only.leaf = FALSE, use.alias = TRUE)
    expect_true(all(grepl("alias",names(actual))))
    
    # signalNode
    actual <- signalNode(tree = tinyTree, node = c("Node_12","Node_15"))
    expect_equal(unname(actual), 12L)
    actual <- signalNode(tree = tinyTree, node = c("Node_13","Node_15"))
    expect_equal(unname(actual), 12L)
    actual <- signalNode(tree = tinyTree, node = c("Node_14","Node_15"))
    expect_equal(unname(actual), c(14L,15L))
    
    # trackNode
    actual <- showNode(trackNode(tree = tinyTree))
    expect_equal(unname(actual),seq_len(19L))
    expect_true(all(grepl("alias",names(actual))))
    
    # unionLeaf 
    actual <- unionLeaf(tree = tinyTree, node = c("Node_16","Node_19"))
    expect_equal(actual, seq.int(4L, 8L))
    actual <- unionLeaf(tree = tinyTree, node = c("t2","Node_19"))
    expect_equal(actual, c(1L, seq.int(7L, 8L)))
    
    # pruneTree
    actual <- pruneTree(tree = tinyTree, rmLeaf = c(4, 5), mergeSingle = FALSE)
    expect_equal(countNode(actual),16L)
    actual <- pruneTree(tree = tinyTree, rmLeaf = c(4, 5), mergeSingle = TRUE)
    expect_equal(countNode(actual),15L)
})
