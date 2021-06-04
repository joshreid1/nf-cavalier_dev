FROM jemunro/igv-snapshot-automator:latest

LABEL \
  author="Jacob Munro" \
  description="Container for cavalier" \
  maintainer="Bahlo Lab"

# Install procps so that Nextflow can poll CPU usage and
# Install recent version of r from cran repo, and r package dependencies
# deep clean the apt cache to reduce image/layer size (from nfcore/base)
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        systemd \
        software-properties-common \
        apt-transport-https \
        ca-certificates \
        procps \
        devscripts \
        libcurl4-openssl-dev \
        build-essential \
        zlib1g-dev \
        libxml2-dev \
        libssl-dev \
        gfortran \
        libblas-dev \
        liblapack-dev \
        libpng-dev \
        libcairo2-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
    && apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/debian buster-cran35/' \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        r-base \
        r-base-core \
        r-recommended \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Install required R packages
COPY inst/install_packages.R  inst/cran_packages.txt inst/bioc_packages.txt /
RUN Rscript --vanilla install_packages.R CRAN:cran_packages.txt BIOC:bioc_packages.txt

# Install cavalier R package
COPY inst/github_packages.txt /
RUN Rscript --vanilla install_packages.R GITHUB:github_packages.txt

# Instruct R processes to use these empty files instead of clashing with a local version (from nfcore/base)
RUN touch .Rprofile .Renviron
ENV R_LIBS_USER=/usr/local/lib/R/site-library TZ=Etc/UTC