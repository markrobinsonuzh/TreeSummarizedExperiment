context("asPhylo")

test_that("asPhylo works properly", {
    df1 <- data.frame(L1 = rep("A", 10),
                     L2 = rep(c("B1", "B2"), c(5, 5)),
                     L3 = rep(NA, 10),
                     L4 = rep(NA, 10),
                     L5 = letters[1:10])
    tr1 <- asPhylo(df1)
    expect_equal(class(tr1), "phylo")
    expect_equal(tr1$Nnode, 3)
    expect_equal(tr1$tip.label, letters[1:10])
    
    # test usNA
    df2 <- data.frame(L1 = rep("A", 10),
                      L2 = rep(c("B1", "B2"), c(5, 5)),
                      L3 = rep("unknown", 10),
                      L4 = rep(c("missing", ""), 5),
                      L5 = letters[1:10])
    tr2 <- asPhylo(df2, asNA = c("unknown", "missing", ""))
    expect_equal(tr1, tr2)    

    df3 <- data.frame(L1 = rep("A", 10),
                      L2 = rep(c("B1", "B2"), c(5, 5)),
                      L5 = letters[1:10])
    tr3 <- asPhylo(df3)
    expect_equal(tr1, tr3)  
    
    # test column_order 
    df4 <- data.frame(L1 = rep("A", 10),
                      L5 = letters[1:10],
                      L2 = rep(c("B1", "B2"), c(5, 5)))
    tr4 <- asPhylo(df3, column_order = c("L1", "L2", "L5"))
    expect_equal(tr3, tr4) 
})