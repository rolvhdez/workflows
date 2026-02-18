# Author: Roberto Olvera Hernandez
# Date: 2026-02-18
FROM rocker/verse:4.3.1

# Install system dependencies (e.g., for 'sf', 'curl', 'xml2' packages) ---
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfreetype6-dev \
    libxml2-dev && \
    rm -rf /var/lib/apt/lists/*

# Copy the list of R packages to install first (for better caching) ---
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile

# Install renv and restore the project library ---
RUN Rscript -e 'renv::restore(lockfile = "renv.lock", repos = NULL)'

# Copy the rest of the pipeline code ---
COPY utils/ utils/
COPY heritability-comparisons.R heritability-comparisons.R
COPY heritability-plots.R heritability-plots.R