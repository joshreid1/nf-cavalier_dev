FROM bahlolab/r-verse-bioc:dev

LABEL \
  author="Jacob Munro" \
  description="Container for cavalier" \
  maintainer="Bahlo Lab"

# Install the R packages
COPY install.R /
RUN Rscript install.R --vanilla
