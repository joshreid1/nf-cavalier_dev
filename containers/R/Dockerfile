FROM rocker/tidyverse:4.0.5

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
        libbz2-dev \
        liblzma-dev
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Install IGV
RUN wget http://data.broadinstitute.org/igv/projects/downloads/2.11/IGV_2.11.0.zip -O IGV.zip \
    && unzip IGV.zip \
    && rm IGV.zip

# Install required R packages
COPY install_packages.R  cran_packages.txt bioc_packages.txt /
RUN Rscript --vanilla install_packages.R CRAN:cran_packages.txt BIOC:bioc_packages.txt

# Install cavalier R package
COPY github_packages.txt /
RUN Rscript --vanilla install_packages.R GITHUB:github_packages.txt

# set R ENV variables and add IGV to R PATH
RUN sed 's:^PATH=:PATH=/IGV_2.11.0\::' -i /usr/local/lib/R/etc/Renviron
ENV PATH="/IGV_2.11.0:${PATH}" \
    TZ=Etc/UTC \
    R_HOME=/usr/local/lib/R/ \
    R_ENVIRON=/usr/local/lib/R/etc/Renviron R_LIBS_USER=/usr/local/lib/R/site-library