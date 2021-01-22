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
    
    
    # findSibling
    expect_equivalent(findSibling(tree = tinyTree, node = c("Node_12","Node_15")),
                      c(10, 13))
    expect_equivalent(findSibling(tree = tinyTree, node = c(10, 18)),
                      c(12,6))
    expect_error(findSibling(tree = tinyTree, node = 11))
    
    expect_equivalent(findSibling(tree = tinyTree, node = 17), 19)
    expect_equivalent(findSibling(tree = tinyTree, node = 16), 9)
    expect_equivalent(findSibling(tree = tinyTree, node = 16:17), c(9, 19))
    expect_error(findSibling(tree = tinyTree, node = 51))
    
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
    
    # joinNode
    actual <- joinNode(tree = tinyTree, node = c("Node_12","Node_15"))
    expect_equal(unname(actual), 12L)
    actual <- joinNode(tree = tinyTree, node = c("Node_13","Node_15"))
    expect_equal(unname(actual), 12L)
    actual <- joinNode(tree = tinyTree, node = c("Node_14","Node_15"))
    expect_equal(unname(actual), c(14L,15L))
    actual <- joinNode(tree = tinyTree, node = c('t4','t9'))
    expect_equal(actual, c("Node_18" = 18))
    actual <- joinNode(tree = tinyTree, node = c('t10','Node_18', 't8'))
    expect_equal(actual, c("Node_17" = 17, "t10" = 7))
    actual <- joinNode(tree = tinyTree, node = c(5, 10))
    expect_equal(actual, c("t4" = 5, "t3" = 10))
    expect_error(joinNode(node = c('t4','t9', 't11'), tree = tinyTree,
                          use.alias = FALSE))
    expect_warning(signalNode(node = 4:5, tree = tinyTree),
                   "'signalNode' is deprecated.")
    # trackNode
    actual <- showNode(trackNode(tree = tinyTree))
    expect_equal(unname(actual),seq_len(19L))
    expect_true(all(grepl("alias",names(actual))))
    
    # unionLeaf 
    actual <- unionLeaf(tree = tinyTree, node = c("Node_16","Node_19"))
    expect_equal(actual, seq.int(4L, 8L))
    actual <- unionLeaf(tree = tinyTree, node = c("t2","Node_19"))
    expect_equal(actual, c(1L, seq.int(7L, 8L)))
    
    # shareNode
    expect_equal(shareNode(tree = tinyTree, node = c('t4','t9')),
                      c("Node_18" = 18))
    expect_equal(shareNode(tree = tinyTree, node = c('t4','t9')),
                      c("Node_18" = 18))
    expect_equal(shareNode(tree = tinyTree, node = c('t4','t9', 't1')),
                 c("Node_16" = 16))
    expect_error(shareNode(node = c('t4','t9', 't11'), tree = tinyTree,
                           use.alias = FALSE))
    
    
    # countNode
    actual <- countNode(tree = tinyTree)
    expect_equal(actual, 19L)
    
    # isLeaf
    actual <- isLeaf(tree = tinyTree, node = 1:10)
    expect_true(all(actual))
    actual <- isLeaf(tree = tinyTree, node = 11)
    expect_false(actual)
    
    # asLeaf
    newTree <- asLeaf(tree = tinyTree, node = "Node_16")
    expect_true(isLeaf(tree = newTree, node = "Node_16"))
    expect_equal(asLeaf(tree = tinyTree, node = "Node_16"),
                asLeaf(tree = tinyTree, node = 16))
    })
