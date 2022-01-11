#'@title DEG_list_extract
#'
#'@description To extract DEG list for further GO, KEGG enrichment analysis and vocanol plot painting
#'
#'@details To extract DEG list for further GO, KEGG enrichment analysis and vocanol plot painting
#'
#'@param DEGs_list needs analyzed list for extraction of DEGs
#'
#'@param log_fold_change_cut_off logFC cut off value. the default value is absolute value 2
#'
#'@param q_value_cut_off FDR cut off value,  the default value is absolute value 0.05
#'
#'@return a list containing 2 lists
#'
#'@importFrom dplyr select arrange %>% contains tibble distinct filter
#'
#'@export
##

DEG_list_extract <- function(DEGs_list, log_fold_change_cut_off = 2, q_value_cut_off = 0.05){
  df_up <- DEGs_list %>%
    filter(adj.P.Val < q_value_cut_off, logFC > log_fold_change_cut_off)
  ids_up <- rownames(df_up)
  df_up <- cbind(df_up, UP = ids_up) %>%
    select(UP)
  rownames(df_up) <- NULL
  df_up <- merge(df_up, MSU, by.x = 'UP', by.y = 'ID') %>%
    select(MSU) %>%
    distinct(MSU) %>%
    filter(MSU != 'None')

  df_down <- DEGs_list %>%
    filter(adj.P.Val < q_value_cut_off, logFC < -log_fold_change_cut_off)
  ids_down <- rownames(df_down)
  df_down <- cbind(df_down, DOWN = ids_down) %>%
    select(DOWN)
  rownames(df_down) <- NULL
  df_down <- merge(df_down, MSU, by.x = 'DOWN', by.y = 'ID') %>%
    select(MSU) %>%
    distinct(MSU) %>%
    filter(MSU != 'None')

    DEG_results <- list(df_up, df_down)
    names(DEG_results) <- c('UP', 'DOWN')

  return(DEG_results)
}



