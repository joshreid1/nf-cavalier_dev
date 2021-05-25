
cran_pkgs <- c('tidyverse', 'matrixStats', 'BiocManager', 'devtools', 'shiny')

install.packages(cran_pkgs,
                 repos='https://cloud.r-project.org',
                 clean = TRUE,
                 verbose = FALSE,
                 quiet = TRUE)

BiocManager::install('org.Hs.eg.db')

devtools::install_github('jemunro/cavalier@latest', force = TRUE)