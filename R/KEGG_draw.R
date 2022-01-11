#'@title KEGG_draw
#'
#'@description Draw KEGG enrichment analysis results
#'
#'@details Draw KEGG enrichment analysis results
#'
#'@param KEGG_list needs KEGG enrichment analysis result
#'
#'@param axis_text_size axis_text_size can adjust the size of words in axis
#'
#'@param axis_title_size axis_title_size can adjust the size of words in title
#'
#'@return a list of 2 graphs, GO_up and GO_down
#'
#'@importFrom dplyr select arrange %>% filter slice
#'
#'@importFrom ggplot2 ggplot geom_bar scale_fill_manual aes theme_bw theme coord_flip labs
#'
#'@export
##

KEGG_draw <- function(KEGG_list, axis_text_size = 14, axis_title_size = 16){
  KEGG_up <- KEGG_list[[1]]
  ids <- rownames(KEGG_up)
  KEGG_up <- cbind(KEGG_up, ids)
  rownames(KEGG_up) <- NULL
  KEGG_up <- KEGG_up %>%
    arrange(desc(gene.number))

  order_KEGG_up <- factor(as.integer(rownames(KEGG_up)), labels = KEGG_up$ids)

  Kdraw_graphic_up <- ggplot(KEGG_up, aes(x = order_KEGG_up, y = gene.number)) +
    geom_bar(stat="identity", width = 0.3, fill = '#fcbd00') +
    theme_bw() +
    theme(text = element_text('sans'),
          axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),legend.position = 'none') +
    coord_flip() +
    labs(x = 'KEGG terms', y = 'Num of Genes', title = 'Up-regulated Genes KEGG enrichment analysis')

  #down
  KEGG_down <- KEGG_list[[2]]
  ids <- rownames(KEGG_down)
  KEGG_down <- cbind(KEGG_down, ids)
  rownames(KEGG_down) <- NULL
  KEGG_down <- KEGG_down %>%
    arrange(desc(gene.number))

  order_KEGG_down <- factor(as.integer(rownames(KEGG_down)), labels = KEGG_down$ids)

  Kdraw_graphic_down <- ggplot(KEGG_down, aes(x = order_KEGG_down, y = gene.number)) +
    geom_bar(stat="identity", width = 0.3, fill = '#468bc9') +
    theme_bw() +
    theme(text = element_text('sans'),
          axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),legend.position = 'none') +
    coord_flip() +
    labs(x = 'KEGG terms', y = 'Num of Genes', title = 'Down-regulated Genes KEGG enrichment analysis')

  graph_list <- list(Kdraw_graphic_up, Kdraw_graphic_down)
  return(graph_list)

}
