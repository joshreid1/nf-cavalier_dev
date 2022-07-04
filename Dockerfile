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
ARG NAME='cavalier_dev'
COPY environment.yml /
RUN conda env create -f /environment.yml \
    && conda clean -a -y \
    && conda env export --name $NAME > $NAME.yml

# Install Cavalier R package
RUN /opt/conda/envs/$NAME/bin/R --slave --vanilla -e \
    "devtools::install_github('jemunro/cavalier@363588d11d982b2909c100dbcd3f06c46ce576d4', \
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