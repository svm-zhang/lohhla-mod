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

RUN R -e "install.packages(c('remotes', 'BiocManager'), repos='http://cran.rstudio.com/')"

RUN mkdir -p /opt/lohhlamod

WORKDIR /opt/lohhlamod
# This way avoids to rebuild dependencies using layer caching
COPY DESCRIPTION .
RUN R -e "remotes::install_deps(dependencies = TRUE, repos = BiocManager::repositories())"

# Only build the library every time changes are made
COPY . .
RUN R -e "remotes::install_local()"

# Call as binaries b/c of shebang
RUN chmod +x /opt/lohhlamod/lohhlamod.R /opt/lohhlamod/lohhlaplot.R
RUN ln -s /opt/lohhlamod/lohhlamod.R /usr/local/bin/lohhlamod
RUN ln -s /opt/lohhlamod/lohhlaplot.R /usr/local/bin/lohhlaplot

# This should match the map point in compose config
WORKDIR /lohhla_runs

