FROM nfcore/base:2.1

LABEL \
  author="Jacob Munro" \
  description="Container for cavalier" \
  maintainer="Bahlo Lab"

# Install xvfb and xauth
# deep clean the apt cache to reduce image/layer size (from nfcore/base)
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        xvfb \
        xauth \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Install the conda environment
ARG NAME='cavalier'
COPY environment.yml /
RUN conda update -n base conda -y \
    && conda install mamba -n base -c conda-forge -y \
    && mamba env create -f /environment.yml \
    && conda clean -a -y \
    && conda env export --name $NAME > $NAME.yml

# Install Cavalier R package
RUN /opt/conda/envs/$NAME/bin/R --slave --vanilla -e \
    "devtools::install_github('jemunro/cavalier@a0f2fd9bbf60aaa11f27734fe131c6fbfc2257e1', \
        force = TRUE, upgrade = 'never')"

# ensure igv.sh exists and add conda executables to R PATH
RUN cp /opt/conda/envs/$NAME/bin/igv /opt/conda/envs/$NAME/bin/igv.sh \
    && echo "PATH=/opt/conda/envs/$NAME/bin:/opt/conda/bin:${PATH}" >> \
    /opt/conda/envs/$NAME/lib/R/etc/Renviron

ENV PATH="/opt/conda/envs/$NAME/bin:/opt/conda/bin:${PATH}" \
    TZ=Etc/UTC \
    R_HOME=/opt/conda/envs/$NAME/lib/R/ \
    R_ENVIRON=/opt/conda/envs/$NAME/lib/R/etc/Renviron \
    R_LIBS_USER=/opt/conda/envs/$NAME/lib/R/site-library