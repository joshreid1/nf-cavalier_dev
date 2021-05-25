FROM stevekm/igv-snapshot-automator:20.11.1

LABEL \
  author="Jacob Munro" \
  description="Container for cavalier" \
  maintainer="Bahlo Lab"

# Install procps so that Nextflow can poll CPU usage and
# Install recent version of r from cran repo
# deep clean the apt cache to reduce image/layer size (from nfcore/base)
RUN apt-get update \
    && apt-get install -y software-properties-common apt-transport-https ca-certificates \
    && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran36/' \
    && apt-get update \
    && apt-get install -y procps devscripts r-base r-base-core r-recommended \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Install the R packages (cavalier)
COPY install.R /
RUN Rscript install.R --vanilla
