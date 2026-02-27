FROM rocker/r-ver:4.1.0

# Install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    samtools \
    python3 \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Install remotes to manage specific versions
RUN R -e "install.packages('remotes', repos='http://cran.rstudio.com/')"

# Install CRAN packages
# Note: we use remotes::install_version for splitstackshape to match your 1.4.8 requirement
RUN R -e "install.packages(c( \
  'argparse', \
  'data.table', \
  'ggplot2', \
  'R.utils', \
  'seqinr', \
  'splitstackshape' \
  ), repos='http://cran.rstudio.com/')"

# Install Bioconductor manager and packages
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')" \
    && R -e "BiocManager::install(c('Biostrings', 'Rsamtools'))"

# Manually create a r_lib folder
RUN mkdir -p /usr/local/lib/R/site-library/lohhlamod

# Manually copy libraries to the r_lib folder
COPY R/lib/*.R /usr/local/lib/R/site-library/lohhlamod/
ENV R_LIBS_USER=/usr/local/lib/R/site-library

# Call as binaries b/c of shebang
COPY R/lohhlamod.R /usr/local/bin/lohhlamod
COPY R/lohhlaplot.R /usr/local/bin/lohhlaplot
RUN chmod +x /usr/local/bin/lohhlamod /usr/local/bin/lohhlaplot

# This should match the map point in compose config
WORKDIR /lohhla_runs

