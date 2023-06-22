FROM mambaorg/micromamba:1.4.4-bionic

LABEL \
  author="Jacob Munro" \
  description="Container for Cavalier" \
  maintainer="Bahlo Lab"

# install os deps
USER root
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    procps \
    xvfb \
    xauth \
  && apt-get clean -y \
  && rm -rf /var/lib/apt/lists/*

# install env with micromamba
COPY environment.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml \ 
  && micromamba clean --all --yes

# Install Cavalier R package
RUN /opt/conda/bin/R --slave --vanilla -e \
    "devtools::install_github('jemunro/cavalier@6e41c2a1838da71cf18c336f47107fb1496dbb23', \
        force = TRUE, upgrade = 'never')"

ENV PATH="/opt/conda/bin:${PATH}" \
    TZ=Etc/UTC \
    R_HOME=/opt/conda/lib/R/ \
    R_ENVIRON=/opt/conda/lib/R/etc/Renviron \
    R_LIBS_USER=/opt/conda/lib/R/site-library