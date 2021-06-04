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
        openjdk-11-jdk-headless \
        openjdk-11-jdk \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Install IGV
RUN wget http://data.broadinstitute.org/igv/projects/downloads/2.9/IGV_2.9.5.zip -O IGV.zip \
    && unzip IGV.zip && rm IGV.zip

# Install required R packages
COPY inst/install_packages.R  inst/cran_packages.txt inst/bioc_packages.txt /
RUN Rscript --vanilla install_packages.R CRAN:cran_packages.txt BIOC:bioc_packages.txt

# Install cavalier R package
COPY inst/github_packages.txt /
RUN Rscript --vanilla install_packages.R GITHUB:github_packages.txt

# Instruct R processes to use these empty files instead of clashing with a local version (from nfcore/base)
RUN touch .Rprofile .Renviron
ENV PATH="/IGV_2.9.5/:${PATH}" R_LIBS_USER=/usr/local/lib/R/site-library TZ=Etc/UTC