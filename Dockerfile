FROM bahlolab/r-verse-bioc:dev

LABEL \
  author="Jacob Munro" \
  description="Container for cavalier" \
  maintainer="Bahlo Lab"

# install IGV-snapshot-automater deps
RUN apt-get update \
    && apt-get install -y wget \
    unzip \
    default-jdk \
    xvfb \
    xorg \
    python \
    make \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# install IGV-snapshot-automater
RUN git clone https://github.com/stevekm/IGV-snapshot-automator.git /IGV-snapshot-automator \
    && cd /IGV-snapshot-automator \
    && make install \
    && printf 'new\ngenome hg19\nexit\n' > /genome.bat \
    && xvfb-run --auto-servernum --server-num=1 igv.sh -b /genome.bat

# Install the R packages (cavalier)
COPY install.R /
RUN Rscript install.R --vanilla
