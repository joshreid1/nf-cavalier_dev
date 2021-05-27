
cran_pkgs <- c('tidyverse', 'matrixStats', 'BiocManager', 'devtools', 'shiny')

install.packages(cran_pkgs,
                 repos='https://cloud.r-project.org',
                 clean = TRUE,
                 verbose = FALSE,
                 quiet = TRUE)

BiocManager::install('org.Hs.eg.db')

devtools::install_github('jemunro/cavalier@latest',
                         ref = '1b6d8f0b52c7b045383aa2d93899ef1b6603ba2c',
                         force = TRUE)