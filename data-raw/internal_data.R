## code to prepare `internal_data` dataset goes here

GPL6864_anno <- read.csv('GPL6864.csv')
GPL8852_anno <- read.csv('GPL8852.csv')
GPL2025_anno <- read.csv('GPL2025.csv')

usethis::use_data(GPL6864_anno, GPL8852_anno, GPL2025_anno, internal = TRUE, overwrite = TRUE)
