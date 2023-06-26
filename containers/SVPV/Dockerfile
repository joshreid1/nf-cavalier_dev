FROM mambaorg/micromamba:1.4.4-bionic

LABEL \
  author="Jacob Munro" \
  description="Container for cavalier SVPV" \
  maintainer="Bahlo Lab"

# install os deps
USER root
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    procps \
    git \
    && apt-get clean -y \
  && rm -rf /var/lib/apt/lists/*

# install env with micromamba
COPY environment.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml \ 
  && micromamba clean --all --yes

# Install SVPV
RUN git clone https://github.com/VCCRI/SVPV.git \
    && echo "python /SVPV/SVPV \$@" > /opt/conda/bin/SVPV \
    && chmod +x /opt/conda/bin/SVPV

ENV PATH="/opt/conda/bin:${PATH}" \
    TZ=Etc/UTC \
    R_HOME=/opt/conda/envs/$NAME/lib/R/ \
    R_ENVIRON=/opt/conda/envs/$NAME/lib/R/etc/Renviron \
    R_LIBS_USER=/opt/conda/envs/$NAME/lib/R/site-library