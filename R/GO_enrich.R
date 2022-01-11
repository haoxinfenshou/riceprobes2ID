#'@title GO_enrich
#'
#'@description Perform GO enrichment analysis
#'
#'@details GO enrichment analysis using fisher's exact test and multi test, respectively
#'
#'@param DEG_results needs analyzed list for extraction of DEGs
#'
#'@return a dataframe
#'
#'@importFrom dplyr select arrange %>% contains distinct filter rename
#'
#'@importFrom tibble add_column
#'
#'@export
##

GO_enrich <- function(DEG_results){
  GO_count <- nrow(GO_term %>%
                     select(Gene) %>%
                     distinct(Gene))
  GO_match <- merge(GO_term, DEG_results[[1]], by.x = 'Gene', by.y = 'MSU')
  GO_num <- nrow(GO_match %>%
                   distinct(Gene))
  GO_number <- as.data.frame(table(GO_match$Gene.set.name))
  genomic_GO_number <- as.data.frame(table(GO_term$Gene.set.name))
  merging1 <- merge(GO_number, genomic_GO_number, by.x = 'Var1', by.y = 'Var1')
  GO_df <- merging1 %>%
    add_column(column_new = GO_num) %>%
    add_column(column_new1 = GO_count)

  rownames(GO_df) <- GO_df$Var1
  GO_df <- GO_df %>%
    select(-Var1) %>%
    add_column(A = GO_df$Freq.x) %>%
    add_column(B = GO_df$Freq.y - GO_df$Freq.x) %>%
    add_column(C = GO_df$column_new - GO_df$Freq.x)

  GO_df <- GO_df %>%
    add_column(D = GO_df$column_new1 - GO_df$A - GO_df$B - GO_df$C) %>%
    select(A, B, C, D)

  GO_df <- GO_df %>%
    add_column(AB = GO_df$A + GO_df$B) %>%
    add_column(CD = GO_df$C + GO_df$D) %>%
    add_column(AC = GO_df$A + GO_df$C) %>%
    add_column(N = 45432)


  Fisher_Exact_Test <- function(row){
    test.mat <- rbind(c(row['A'], row['B']),
                      c(row['C'], row['D']))
    test.results <- fisher.test(test.mat, alternative = "greater")
    return(test.results$p.value)
  }

  fisher.exact.test <- apply(X = GO_df, MARGIN = 1, FUN = Fisher_Exact_Test)

  fisher.exact.test.q <- p.adjust(p = fisher.exact.test,
                                  method = 'BH')
  out_up_GO <- cbind(GO_df, fisher.exact.test, fisher.exact.test.q) %>%
    select(A, AC, AB, N, fisher.exact.test, fisher.exact.test.q) %>%
    rename(gene.number = A, gene.list = AC, gene.number.in.GO = AB, total.number = N,
           p.value = fisher.exact.test, FDR = fisher.exact.test.q) %>%
    filter(FDR < 0.05)

  # down
  GO_count <- nrow(GO_term %>%
                     select(Gene) %>%
                     distinct(Gene))
  GO_match <- merge(GO_term, DEG_results[[2]], by.x = 'Gene', by.y = 'MSU')
  GO_num <- nrow(GO_match %>%
                   distinct(Gene))
  GO_number <- as.data.frame(table(GO_match$Gene.set.name))
  genomic_GO_number <- as.data.frame(table(GO_term$Gene.set.name))
  merging1 <- merge(GO_number, genomic_GO_number, by.x = 'Var1', by.y = 'Var1')
  GO_df <- merging1 %>%
    add_column(column_new = GO_num) %>%
    add_column(column_new1 = GO_count)

  rownames(GO_df) <- GO_df$Var1
  GO_df <- GO_df %>%
    select(-Var1) %>%
    add_column(A = GO_df$Freq.x) %>%
    add_column(B = GO_df$Freq.y - GO_df$Freq.x) %>%
    add_column(C = GO_df$column_new - GO_df$Freq.x)

  GO_df <- GO_df %>%
    add_column(D = GO_df$column_new1 - GO_df$A - GO_df$B - GO_df$C) %>%
    select(A, B, C, D)

  GO_df <- GO_df %>%
    add_column(AB = GO_df$A + GO_df$B) %>%
    add_column(CD = GO_df$C + GO_df$D) %>%
    add_column(AC = GO_df$A + GO_df$C) %>%
    add_column(N = 45432)


  Fisher_Exact_Test <- function(row){
    test.mat <- rbind(c(row['A'], row['B']),
                      c(row['C'], row['D']))
    test.results <- fisher.test(test.mat, alternative = "greater")
    return(test.results$p.value)
  }

  fisher.exact.test <- apply(X = GO_df, MARGIN = 1, FUN = Fisher_Exact_Test)

  fisher.exact.test.q <- p.adjust(p = fisher.exact.test,
                                  method = 'BH')
  out_down_GO <- cbind(GO_df, fisher.exact.test, fisher.exact.test.q) %>%
    select(A, AC, AB, N, fisher.exact.test, fisher.exact.test.q) %>%
    rename(gene.number = A, gene.list = AC, gene.number.in.GO = AB, total.number = N,
           p.value = fisher.exact.test, FDR = fisher.exact.test.q) %>%
    filter(FDR < 0.05)

  GO_result <- list(out_up_GO, out_down_GO)
  names(GO_result) <- c('UP_GO', 'DOWN_GO')
  return(GO_result)

}
