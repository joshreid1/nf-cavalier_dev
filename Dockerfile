FROM rocker/tidyverse:3.6.3

LABEL \
  author="Jacob Munro" \
  description="Container for cavalier" \
  maintainer="Bahlo Lab"

# Install procps so that Nextflow can poll CPU usage and
# Install recent version of r from cran repo, and r package dependencies
# deep clean the apt cache to reduce image/layer size (from nfcore/base)
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        xvfb \
        xauth \
        openjdk-11-jdk-headless \
        openjdk-11-jdk \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Install IGV
RUN wget http://data.broadinstitute.org/igv/projects/downloads/2.9/IGV_2.9.5.zip -O IGV.zip \
    && unzip IGV.zip \
    && rm IGV.zip
# # and download genomes
#    && mkdir /igv \
#    && printf '%s\n%s\n' '#!/bin/bash' '/IGV_2.9.5/igv.sh --igvDirectory /igv "$@"' > /igv/igv.sh \
#    && chmod +x /igv/igv.sh \
#    && printf '%s\n%s\n%s\n%s\n' 'new' 'genome hg19' 'genome hg38' 'exit' > /genome.bat \
#    && xvfb-run --auto-servernum --server-num=1 /igv/igv.sh -b /genome.bat

# Install required R packages
COPY inst/install_packages.R  inst/cran_packages.txt inst/bioc_packages.txt /
RUN Rscript --vanilla install_packages.R CRAN:cran_packages.txt BIOC:bioc_packages.txt

# Install cavalier R package
COPY inst/github_packages.txt /
RUN Rscript --vanilla install_packages.R GITHUB:github_packages.txt

# set R ENV variables and add IGV to R PATH
RUN sed 's:^PATH=:PATH=/IGV_2.9.5\::' -i /usr/local/lib/R/etc/Renviron
ENV PATH="/IGV_2.9.5:${PATH}" TZ=Etc/UTC \
    R_HOME=/usr/local/lib/R/ R_ENVIRON=/usr/local/lib/R/etc/Renviron R_LIBS_USER=/usr/local/lib/R/site-library