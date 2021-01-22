context("toTree")

df <- data.frame(R1 = rep("A", 9),
                 R2 = c("B1", rep(c("B2", ""), each =4)),
                 R3 = c("C1", "C2", "C3", "C4", "C5", rep("", 4)),
                 R4 = c(rep(NA, 5), paste0("D", 1:4)))

test_that("toTree could convert a data frame to a phylo", {
    df_complete <- df[1:5, 1:3]
    expect_phylo_output <- function(df, truth) {
        out <- toTree(data = df)
        expect_setequal(class(out), truth)
        lab <- showNode(tree = out, only.leaf = FALSE)
        out <- addLabel(tree = out, label = as.character(lab), on = "all")
        expect_equal(out$tip.label,as.character(seq_len(length(out))))
        expect_equal(out,addLabel(out,as.character(seq_len(length(out) + length(out$node.label)))))
    }
    expect_phylo_nodes <- function(df, truth, use_leaf) {
        out <- toTree(data = df)
        if (use_leaf) {
            nod <- out$tip.label
        } else {
            nod <- out$node.label
        }
        
        expect_setequal(nod, truth)
    }
    
    expect_phylo_output(df_complete, "phylo")
    expect_phylo_nodes(df_complete, 
                       c("C1", "C2", "C3", "C4", "C5"), 
                       use_leaf = TRUE)
    expect_phylo_nodes(df_complete, 
                       c("R1:A",  "R2:B1", "R2:B2"), 
                       use_leaf = FALSE)
    
    expect_phylo_nodes(df[5:9, ], 
                       c("R1:A", "R2:", "R3:", "R2:B2", "R3:C5"), 
                       use_leaf = FALSE)
    expect_phylo_nodes(df[5:9, ], 
                       c("NA", "D1", "D2", "D3", "D4"), 
                       use_leaf = TRUE)
})

test_that("toTree give warnings", {
    
    expect_error(toTree(data = df))
    expect_warning(toTree(data = df[, 1:2]))
    
})


test_that("toTree loop resolving", {
    actual <- expect_warning(detectLoop(tax_tab = df),
                             "Detected R4 from different R3")
    expect_s3_class(actual,"data.frame")
    expect_equal(nrow(actual),5L)
    expect_true(all(is.na(actual$child)))
    
    df <- data.frame(A = rep("a", 8),
                     B = rep (c("b1", "b2", "b3", "b4"), each = 2),
                     C = paste0("c", c(1, 2, 2, 3:7)),
                     D = paste0("d", 1:8))
    actual <- expect_warning(detectLoop(tax_tab = df),
                             "Detected C from different B")
    expect_s3_class(actual,"data.frame")
    expect_equal(nrow(actual),2L)
    expect_true(all(actual$child == "c2"))
    
    actual <- expect_warning(resolveLoop(tax_tab = df),
                             "Detected C from different B")
    expect_equal(nrow(df),nrow(actual))
    expect_equal(actual$C[2L],"c2_1")
    expect_equal(actual$C[3L],"c2_2")
})
