
cran_pkgs <- c('tidyverse', 'matrixStats', 'BiocManager')

install.packages(cran_pkgs, repos='https://cloud.r-project.org')

BiocManager::install('org.Hs.eg.db')

devtools::install_github('jemunro/cavalier@latest', force = TRUE)