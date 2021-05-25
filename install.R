
cran_pkgs <- c('cowplot', 'ggforce', 'matrixStats', 'BiocManager')

install.packages(cran_pkgs, repos='https://cloud.r-project.org')

BiocManager::install('org.Hs.eg.db')

devtools::install_github('jemunro/cavalier@latest', force = TRUE)