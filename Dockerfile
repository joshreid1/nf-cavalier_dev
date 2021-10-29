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
COPY inst/environment.yml /
RUN conda env create -f /environment.yml \
    && conda clean -a -y \
    && conda env export --name $NAME > $NAME.yml

# Install cavalier R package
COPY inst/install_packages.R inst/github_packages.txt /
RUN /opt/conda/envs/$NAME/bin/Rscript --vanilla install_packages.R GITHUB:github_packages.txt

# set R ENV variables and add IGV to R PATH
RUN cp /opt/conda/envs/$NAME/bin/igv /opt/conda/envs/$NAME/bin/igv.sh \
    && echo 'PATH=/opt/conda/envs/$NAME/bin' >> /opt/conda/envs/$NAME/lib/R/etc/Renviron
ENV PATH="/opt/conda/envs/$NAME/bin:/opt/conda/bin:${PATH}" \
    TZ=Etc/UTC \
    R_HOME=/opt/conda/envs/$NAME/lib/R/ \
    R_ENVIRON=/opt/conda/envs/$NAME/lib/R/etc/Renviron \
    R_LIBS_USER=/opt/conda/envs/$NAME/lib/R/site-library