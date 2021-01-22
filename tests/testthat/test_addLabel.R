context("addLabel")

set.seed(1)
n <- 10
tr <- ape::rtree(n)

test_that("addLabel works properly", {
    tr1 <- addLabel(tree = tr, on = "leaf")
    tr2 <- addLabel(tree = tr, on = "internal")
    tr3 <- addLabel(tree = tr, on = "all")
    
    expect_equal(tr1$tip.label, paste0("Node_", seq_len(n)))
    expect_equal(tr1$tip.label, tr3$tip.label)
    expect_equal(tr2$node.label, tr3$node.label)
    expect_false(setequal(tr1$tip.label, tr2$tip.label))
    
    tr4 <- addLabel(tree = tr, on = "leaf", label = LETTERS[seq_len(n)])
    tr5 <- addLabel(tree = tr, on = "internal", label = LETTERS[seq_len(n-1)])
    tr6 <- addLabel(tree = tr, on = "all", label = LETTERS[seq_len(2*n-1)])
    expect_equal(tr4$tip.label, LETTERS[seq_len(n)])
    expect_equal(tr5$node.label, LETTERS[seq_len(n-1)])
    expect_equal(c(tr6$tip.label, tr6$node.label), LETTERS[seq_len(2*n-1)])
})

