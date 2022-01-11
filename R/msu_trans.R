#'@title msu_trans
#'
#'@description Convert RAP-ID to MSU-ID
#'
#'@details Convert RAP-ID to MSU-ID
#'
#'@param DEGs_list needs a method to DEGs_list and the rownames should be RAP-ID for MSU ID conversion
#'
#'@return a dataframe
#'
#'@importFrom dplyr select %>%
#'
#'@export
##

msu_trans <- function(DEGs_list){
  ID <- rownames(DEGs_list)
  bind1 <- cbind(ID, DEGs_list)
  bind2 <- merge(bind1, MSU, by.x = 'ID', by.y = 'ID')
  row.names(bind2) <- bind2[, 1]
  bind3 <- bind2 %>%
    select(-ID)
  return(bind3)
}
