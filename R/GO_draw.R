#'@title GO_draw
#'
#'@description Draw GO enrichment analysis results
#'
#'@details Draw GO enrichment analysis results
#'
#'@param GO_list needs GO enrichment analysis result
#'
#'@param term_num needs a cut off number. To draw a tidy graph, only top 15 terms in MF, CC and BP, respectively, will be selected to draw the graph. 15 is default number which users can set the num.
#'
#'@param gene_num minimum number of gene in one GO term, the default value is 3
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

GO_draw <- function(GO_list, term_num = 15, gene_num = 3, axis_text_size = 12, axis_title_size = 14){
  GO_up <- GO_list[[1]]
  ids <- rownames(GO_up)
  GO_up <- cbind(GO_up, ids)
  GO_up <- merge(GO_up, list, by.x = 'ids', by.y = 'Gene.set.name') %>%
    select(-X)

  MF_up <- GO_up %>%
    subset(Category == 'MF') %>%
    arrange(desc(gene.number)) %>%
    slice(1:term_num)
  BP_up <- GO_up %>%
    subset(Category == 'BP') %>%
    arrange(desc(gene.number)) %>%
    slice(1:term_num)
  CC_up <- GO_up %>%
    subset(Category == 'CC') %>%
    arrange(desc(gene.number)) %>%
    slice(1:term_num)

  GO_up <- rbind(MF_up, BP_up, CC_up) %>%
    filter(gene.number > gene_num)

  order_GO1 <- factor(as.integer(rownames(GO_up)), labels = GO_up$ids)

  draw_graphic_up <- ggplot(GO_up, aes(x = order_GO1, y = gene.number, fill = Category)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c('#468bc9', '#f84a2b', '#fcbd00')) +
    theme_bw() +
    theme(text = element_text('sans'),
          axis.text = element_text(size = axis_text_size, colour = 'black'),
          axis.title = element_text(size = axis_title_size)) +
    coord_flip() +
    labs(x = 'GO terms', y = 'Num of Genes', title = 'Up-regulated Genes GO enrichment analysis')
  #down
  GO_down <- GO_list[[2]]
  ids <- rownames(GO_down)
  GO_down <- cbind(GO_down, ids)
  GO_down <- merge(GO_down, list, by.x = 'ids', by.y = 'Gene.set.name') %>%
    select(-X)

  MF_down <- GO_down %>%
    subset(Category == 'MF') %>%
    arrange(desc(gene.number)) %>%
    slice(1:term_num)
  BP_down <- GO_down %>%
    subset(Category == 'BP') %>%
    arrange(desc(gene.number)) %>%
    slice(1:term_num)
  CC_down <- GO_down %>%
    subset(Category == 'CC') %>%
    arrange(desc(gene.number)) %>%
    slice(1:term_num)

  GO_down <- rbind(MF_down, BP_down, CC_down) %>%
    filter(gene.number > gene_num)

  order_GO <- factor(as.integer(rownames(GO_down)), labels = GO_down$ids)

  draw_graphic_down <- ggplot(GO_down, aes(x = order_GO, y = gene.number, fill = Category)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c('#468bc9', '#f84a2b', '#fcbd00')) +
    theme_bw() +
    theme(text = element_text('sans'),
          axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size)) +
    coord_flip() +
    labs(x = 'GO terms', y = 'Num of Genes', title = 'Down-regulated Genes GO enrichment analysis')

  graph_list <- list(draw_graphic_up, draw_graphic_down)

}
