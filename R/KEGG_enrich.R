#'@title KEGG_enrich
#'
#'@description Perform KEGG enrichment analysis
#'
#'@details KEGG enrichment analysis using fisher's exact test and multi test, respectively
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

KEGG_enrich <- function(DEG_results){
  KEGG_count <- nrow(KEGG_term %>%
                     select(Gene) %>%
                     distinct(Gene))
  KEGG_match <- merge(KEGG_term, DEG_results[[1]], by.x = 'Gene', by.y = 'MSU')
  KEGG_num <- nrow(KEGG_match %>%
                   distinct(Gene))
  KEGG_number <- as.data.frame(table(KEGG_match$Gene.Set.Name))
  genomic_KEGG_number <- as.data.frame(table(KEGG_term$Gene.Set.Name))
  merging1 <- merge(KEGG_number, genomic_KEGG_number, by.x = 'Var1', by.y = 'Var1')
  KEGG_df <- merging1 %>%
    add_column(column_new = KEGG_num) %>%
    add_column(column_new1 = KEGG_count)

  rownames(KEGG_df) <- KEGG_df$Var1
  KEGG_df <- KEGG_df %>%
    select(-Var1) %>%
    add_column(A = KEGG_df$Freq.x) %>%
    add_column(B = KEGG_df$Freq.y - KEGG_df$Freq.x) %>%
    add_column(C = KEGG_df$column_new - KEGG_df$Freq.x)

  KEGG_df <- KEGG_df %>%
    add_column(D = KEGG_df$column_new1 - KEGG_df$A - KEGG_df$B - KEGG_df$C) %>%
    select(A, B, C, D)

  KEGG_df <- KEGG_df %>%
    add_column(AB = KEGG_df$A + KEGG_df$B) %>%
    add_column(CD = KEGG_df$C + KEGG_df$D) %>%
    add_column(AC = KEGG_df$A + KEGG_df$C) %>%
    add_column(N = 45432)
  Fisher_Exact_Test <- function(row){
    test.mat <- rbind(c(row['A'], row['B']),
                      c(row['C'], row['D']))
    test.results <- fisher.test(test.mat, alternative = "greater")
    return(test.results$p.value)
  }

  fisher.exact.test <- apply(X = KEGG_df, MARGIN = 1, FUN = Fisher_Exact_Test)

  fisher.exact.test.q <- p.adjust(p = fisher.exact.test,
                                  method = 'BH')
  out_up_KEGG <- cbind(KEGG_df, fisher.exact.test, fisher.exact.test.q) %>%
    select(A, AC, AB, N, fisher.exact.test, fisher.exact.test.q) %>%
    rename(gene.number = A, gene.list = AC, gene.number.in.KEGG = AB, total.number = N,
           p.value = fisher.exact.test, FDR = fisher.exact.test.q) %>%
    filter(FDR < 0.05)

  # down
  KEGG_match <- merge(KEGG_term, DEG_results[[2]], by.x = 'Gene', by.y = 'MSU')
  KEGG_num <- nrow(KEGG_match %>%
                     distinct(Gene))
  KEGG_number <- as.data.frame(table(KEGG_match$Gene.Set.Name))
  genomic_KEGG_number <- as.data.frame(table(KEGG_term$Gene.Set.Name))
  merging1 <- merge(KEGG_number, genomic_KEGG_number, by.x = 'Var1', by.y = 'Var1')
  KEGG_df <- merging1 %>%
    add_column(column_new = KEGG_num) %>%
    add_column(column_new1 = KEGG_count)

  rownames(KEGG_df) <- KEGG_df$Var1
  KEGG_df <- KEGG_df %>%
    select(-Var1) %>%
    add_column(A = KEGG_df$Freq.x) %>%
    add_column(B = KEGG_df$Freq.y - KEGG_df$Freq.x) %>%
    add_column(C = KEGG_df$column_new - KEGG_df$Freq.x)

  KEGG_df <- KEGG_df %>%
    add_column(D = KEGG_df$column_new1 - KEGG_df$A - KEGG_df$B - KEGG_df$C) %>%
    select(A, B, C, D)

  KEGG_df <- KEGG_df %>%
    add_column(AB = KEGG_df$A + KEGG_df$B) %>%
    add_column(CD = KEGG_df$C + KEGG_df$D) %>%
    add_column(AC = KEGG_df$A + KEGG_df$C) %>%
    add_column(N = 45432)
  Fisher_Exact_Test <- function(row){
    test.mat <- rbind(c(row['A'], row['B']),
                      c(row['C'], row['D']))
    test.results <- fisher.test(test.mat, alternative = "greater")
    return(test.results$p.value)
  }

  fisher.exact.test <- apply(X = KEGG_df, MARGIN = 1, FUN = Fisher_Exact_Test)

  fisher.exact.test.q <- p.adjust(p = fisher.exact.test,
                                  method = 'BH')
  out_down_KEGG <- cbind(KEGG_df, fisher.exact.test, fisher.exact.test.q) %>%
    select(A, AC, AB, N, fisher.exact.test, fisher.exact.test.q) %>%
    rename(gene.number = A, gene.list = AC, gene.number.in.KEGG = AB, total.number = N,
           p.value = fisher.exact.test, FDR = fisher.exact.test.q) %>%
    filter(FDR < 0.05)

  KEGG_result <- list(out_up_KEGG, out_down_KEGG)
  names(KEGG_result) <- c('UP_KEGG', 'DOWN_KEGG')
  return(KEGG_result)

}
