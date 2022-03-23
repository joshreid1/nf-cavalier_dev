FROM nfcore/base:2.1

LABEL \
  author="Jacob Munro" \
  description="Container for SVPV" \
  maintainer="Bahlo Lab"

# Install the conda environment
ARG NAME='svpv'
COPY environment.yml /
RUN conda update -n base conda -y \
    && conda env create -f /environment.yml \
    && conda clean -a -y \
    && conda env export --name $NAME > $NAME.yml

# Install SVPV
RUN git clone https://github.com/VCCRI/SVPV.git \
    && echo "python /SVPV/SVPV \$@" > /opt/conda/bin/SVPV \
    && chmod +x /opt/conda/bin/SVPV

ENV PATH="/opt/conda/envs/$NAME/bin:/opt/conda/bin:${PATH}" \
    TZ=Etc/UTC \
    R_HOME=/opt/conda/envs/$NAME/lib/R/ \
    R_ENVIRON=/opt/conda/envs/$NAME/lib/R/etc/Renviron \
    R_LIBS_USER=/opt/conda/envs/$NAME/lib/R/site-library