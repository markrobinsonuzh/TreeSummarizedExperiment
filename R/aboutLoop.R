#' Detect loops 
#' \code{detectLoop} detects loops
#' 
#' @param tax_tab a data frame where columns store hierarchical levels. The
#'   columns from the left to the right correspond nodes from the root to the
#'   leaf.
#' @importFrom dplyr '%>%' mutate_if
#' @export
#' @author Ruizhu Huang
#' @return a data frame
#' @examples 
#' 
#' df <- data.frame(A = rep("a", 8),
#'                  B = rep (c("b1", "b2", "b3", "b4"), each = 2),
#'                  C = paste0("c", c(1, 2, 2, 3:7)),
#'                  D = paste0("d", 1:8))
#'  
#' # The result means that a loop is caused by 'b1' and 'b2' in column 'B' and
#' # 'c2' in column 'C' (a-b1-c2; a-b2-c2)
#' detectLoop(tax_tab = df)
#' 
#' df <- data.frame(R1 = rep("A", 6),
#'                      R2 = c("B1", rep("B2", 4), "B3"),
#'                      R3 = c("C1", "C2", "C3", NA, NA, NA),
#'                      R4 = c("D1", "D2", "D3", NA, NA, NA),
#'                      R5 = paste0("E", 1:6))
#' detectLoop(tax_tab = df)
#' 
#' df <- data.frame(R1 = rep("A", 7),
#'                      R2 = c("B1", rep("B2", 4), "B3", "B3"),
#'                      R3 = c("C1", "C2", "C3", "", "", "", ""),
#'                      R4 = c("D1", "D2", "D3", "", "", "", ""),
#'                      R5 = paste0("E", 1:7))
#' detectLoop(tax_tab = df)
#' 
#' df <- data.frame(R1 = rep("A", 7),
#'                      R2 = c("B1", rep("B2", 4), "B3", "B3"),
#'                      R3 = c("C1", "C2", "C3", NA, NA, NA, NA),
#'                      R4 = c("D1", "D2", "D3", NA, NA, NA, NA),
#'                      R5 = paste0("E", 1:7))
#' detectLoop(tax_tab = df) 

detectLoop <- function(tax_tab){
    
    tax_tab <- tax_tab %>%
        mutate_if(is.factor, as.character) 

    nc <- ncol(tax_tab)
    cnam <- colnames(tax_tab)
    
    tab_list <- lapply(seq_len(nc-1), function(x){
        xx <- table(tax_tab[, x], tax_tab[, x+1], useNA = "always")
        sx <- colSums(xx > 0)
        nx <- names(sx)[sx > 1]
        if (length(nx)) {
            warning("Detected ", cnam[x+1], " from different ",
                    cnam[x], " : ", paste(nx, collapse = " "))
            
            ssx <- xx[, sx>1, drop = FALSE]
            ssx <- ssx[rowSums(ssx) > 0, , drop = FALSE]
            ind <- which(ssx > 0, arr.ind = TRUE)
            
            data.frame(parent = rownames(ssx)[ind[, "row"]],
                       child = colnames(ssx)[ind[, "col"]],
                       parent_column = cnam[x],
                       child_column = cnam[x+1])
        } else {
            NULL
        }
        
    })
    
    out <- do.call(rbind, tab_list)
    data.frame(lapply(out, as.character), stringsAsFactors = FALSE)
}




#' Resolve loops 
#' \code{resolveLoop} resolve loops by adding suffix to the child node. The
#' suffix is "_i" where 'i' is a number. Please see examples.
#' 
#' @param tax_tab a data frame where columns store hierarchical levels. The
#'   columns from the left to the right correspond nodes from the root to the
#'   leaf.
#' @importFrom dplyr '%>%' mutate_if group_by mutate row_number
#' @importFrom rlang .data
#' @export
#' @author Ruizhu Huang
#' @return a data frame
#' @examples 
#' # example 1
#' df <- data.frame(A = rep("a", 8),
#'                  B = rep (c("b1", "b2", "b3", "b4"), each = 2),
#'                  C = paste0("c", c(1, 2, 2, 3:7)),
#'                  D = paste0("d", 1:8))
#'  
#' # The result means that a loop is caused by 'b1' and 'b2' in column 'B' and
#' # 'c2' in column 'C' (a-b1-c2; a-b2-c2)
#' resolveLoop(tax_tab = df)
#' 
#' # example 2
#' taxTab <- data.frame(R1 = rep("A", 5),
#'                      R2 = c("B1", rep("B2", 3), ""),
#'                      R3 = c("C1", "C2", "C3", "", ""),
#'                      R4 = c("D1", "D2", "D3", "", ""),
#'                      R5 = paste0("E", 1:5))
#' 
#' resolveLoop(tax_tab = taxTab)
#' # example 3
#' taxTab <- data.frame(R1 = rep("A", 6),
#'                      R2 = c("B1", rep("B2", 4), ""),
#'                      R3 = c("C1", "C2", "C3", "", "", ""),
#'                      R4 = c("D1", "D2", "D3", "", "", ""),
#'                      R5 = paste0("E", 1:6))
#' 
#' resolveLoop(tax_tab = taxTab)
#' 
#'  # example 3
#' taxTab <- data.frame(
#'             R1 = rep("A", 5),
#'             R2 = c("B1", rep("B2", 3), "B3"),
#'             R3 = c("C1", "C2", "C3", NA, NA),
#'             R4 = c("D1", "D2", "D3", NA, NA),
#'             R5 = paste0("E", 1:5))
#' resolveLoop(tax_tab = taxTab)
#' 
#' 
resolveLoop <- function(tax_tab) {
    # solve the loop by adding suffix
    tax_s <- .solveLoop(tax_tab)
    
    # test whether all loops are solved
    loop_s <- detectLoop(tax_tab = tax_s)
    while (nrow(loop_s)) {
        tax_s <- .solveLoop(tax_s)
        loop_s <- detectLoop(tax_tab = tax_s)
    }
    return(tax_s)
}

.solveLoop <- function(tax_tab){
    # change factor to character
    tax_tab <- tax_tab %>%
        mutate_if(is.factor, as.character) 
    
    df <- detectLoop(tax_tab = tax_tab) 
    if (nrow(df)) {
        df <- df %>% group_by(.data$child) %>% mutate(count = row_number())
        df_list <- split(df, seq(nrow(df)))
        
        tax_x <- tax_tab
        for (i in seq_along(df_list)){
            df_i <- df_list[[i]]
            ind_c <- tax_x[[df_i$child_column]] %in% df_i$child
            ind_p <- tax_x[[df_i$parent_column]] == df_i$parent
            ind <- which(ind_p & ind_c)
            tax_x[ind, df_i$child_column] <- paste(tax_x[ind, df_i$child_column], 
                                                   df_i$count, sep = "_")
        }
        tax_x
    } else {
        tax_tab
    }
    

}
