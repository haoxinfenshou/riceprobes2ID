#'@title volcano_plot
#'
#'@description Draw volcano plot
#'
#'@details Draw volcano plot
#'
#'@param DEG needs analyzed list for drawing volcano plot
#'
#'@param logFC needs cutoff log2 fold change value, the default value is absolute value 2
#'
#'@param q.value needs cutoff q.value, the default value is absolute value 0.05
#'
#'@return a graph
#'
#'@importFrom ggplot2 ggplot aes geom_point scale_color_manual xlim element_blank element_text ylim geom_vline geom_hline theme_bw theme labs
#'
#'@export
##

volcano_plot <- function(DEG, logFC = 2, q.value = 0.05){
  DEG$group <- ifelse(DEG$logFC >= logFC & DEG$adj.P.Val <= q.value, 'UP',
                      ifelse(DEG$logFC <= -logFC & DEG$adj.P.Val <= q.value, 'DOWN', 'NOT SIG'))
  draw_vol <- ggplot(DEG, aes(x = logFC, y = -log10(adj.P.Val))) +
    labs(x = 'log2FC', y = '-log10FDR', title = 'Volcano Plot') +
    geom_point(aes(color = group), size = 3, alpha = 0.5) +
    scale_color_manual(values = c('#468bc9','gray', '#fcbd00')) +
    xlim(c(-max(abs(DEG$logFC)), max(abs(DEG$logFC)))) +
    ylim(c(0, max(-log10(DEG$adj.P.Val)))) +
    geom_vline(xintercept = c(-logFC, logFC), lty=4, lwd = 0.8) +
    geom_hline(yintercept = -log10(q.value), lty=4, lwd = 0.8) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(axis.text = element_text(size = 16, colour = 'black'),
          axis.title = element_text(size = 16), title = element_text(size = 13),
          legend.title = element_blank(),plot.title = element_text(size = 20 ,hjust = 0.5),
          legend.text = element_text(size = 13, colour = 'black')) +
    theme(text = element_text('sans'))

  return(draw_vol)
}
