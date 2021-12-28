#'@title riceprobes2ID
#'
#'@description Convert rice microarray ID to RAP-ID and merge multi probes for one genes
#'
#'@details Convert rice microarray ID to RAP-ID and merge multi probes for one genes
#'
#'@param GPL needs a GPL ID. 2021-12-27, 3 rice microarray platforms are available: GPL6864, GPL8852 and GPL2025
#'
#'@param expr_matrix needs a dataframe, rownames should be microarray ID and colnames should be sample GSM ID
#'
#'@param merge_by needs a method to merge multi probes and now only mean value is available
#'
#'@return a dataframe
#'
#'@importFrom dplyr select arrange %>% contains
#'
#'@export
##

probe_trans <- function(GPL, expr_matrix, merge_by = 'mean'){
  if(GPL == 'GPL6864'){
    ids <- rownames(expr_matrix)
    expr1 <- cbind(expr_matrix, ID = ids)
    merge1 <- merge(expr1, GPL6864_anno, by.x = 'ID', by.y = 'ID')
    clean1 <- merge1 %>%
      select(-ID)
    expr2 <- as.data.frame(lapply(clean1[1:ncol(clean1) - 1], as.numeric))
    expr3 <- cbind(merge1$ACC, expr2)
    expr4 <- expr3 %>%
      select(contains('ACC') | contains('GSM')) %>%
      arrange(merge1$ACC)
    expr5 <- aggregate(x = expr4[,2:ncol(expr4)],
                       by = list(expr4$`merge1$ACC`),
                       FUN = merge_by)
    rownames(expr5) <- expr5[, 1]
    return(expr5[, -1])
  }
  if (GPL == 'GPL8852'){
    ids <- rownames(expr_matrix)
    expr1 <- cbind(expr_matrix, ID = ids)
    merge1 <- merge(expr1, GPL8852_anno, by.x = 'ID', by.y = 'ID')
    clean1 <- merge1 %>%
      select(-ID)
    expr2 <- as.data.frame(lapply(clean1[1:ncol(clean1) - 1], as.numeric))
    expr3 <- cbind(merge1$ACC, expr2)
    expr4 <- expr3 %>%
      select(contains('ACC') | contains('GSM')) %>%
      arrange(merge1$ACC)
    expr5 <- aggregate(x = expr4[,2:ncol(expr4)],
                       by = list(expr4$`merge1$ACC`),
                       FUN = merge_by)
    rownames(expr5) <- expr5[, 1]
    return(expr5[, -1])
  }
  if (GPL == 'GPL2025'){
    ids <- rownames(expr_matrix)
    expr1 <- cbind(expr_matrix, ID = ids)
    merge1 <- merge(expr1, GPL2025_anno, by.x = 'ID', by.y = 'ID')
    clean1 <- merge1 %>%
      select(-ID)
    expr2 <- as.data.frame(lapply(clean1[1:ncol(clean1) - 1], as.numeric))
    expr3 <- cbind(merge1$ACC, expr2)
    expr4 <- expr3 %>%
      select(contains('ACC') | contains('GSM')) %>%
      arrange(merge1$ACC)
    expr5 <- aggregate(x = expr4[,2:ncol(expr4)],
                       by = list(expr4$`merge1$ACC`),
                       FUN = merge_by)
    rownames(expr5) <- expr5[, 1]
    return(expr5[, -1])
  }
}
