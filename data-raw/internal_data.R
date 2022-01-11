## code to prepare `internal_data` dataset goes here

GPL6864_anno <- read.csv('GPL6864.csv')
GPL8852_anno <- read.csv('GPL8852.csv')
GPL2025_anno <- read.csv('GPL2025.csv')
MSU <- read.csv('RAP-MSU.csv')
GO_term <- read.csv('Osa_GO.csv')
KEGG_term <- read.csv('Osa_KEGG.csv')
list <- read.csv('list.csv')

usethis::use_data(GPL6864_anno, GPL8852_anno, GPL2025_anno, MSU, GO_term, KEGG_term, list, internal = TRUE, overwrite = TRUE)

