FROM mambaorg/micromamba:1.5.8-noble

LABEL \
  author="Jacob Munro" \
  description="Container for cavalier VCF utils" \
  maintainer="Bahlo Lab"

# install os deps
USER root
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    procps \
  && apt-get clean -y \
  && rm -rf /var/lib/apt/lists/*

# install env with micromamba
COPY environment.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml \ 
  && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}" \
    TZ=Etc/UTC \
    R_HOME=/opt/conda/lib/R/ \
    R_ENVIRON=/opt/conda/lib/R/etc/Renviron \
    R_LIBS_USER=/opt/conda/lib/R/site-library