# Author: Roberto Olvera Hernandez
# Date: 2025-10-14
FROM rolvhdez/analysis-room:v1.0.0

WORKDIR /opt/

# Install system dependencies (e.g., for 'sf', 'curl', 'xml2' packages) ---
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev && \
    rm -rf /var/lib/apt/lists/*

# Copy the list of R packages to install first (for better caching) ---
# COPY renv/activate.R renv/activate.R
# COPY renv/settings.json renv/settings.json
# COPY renv.lock renv.lock
# COPY .Rprofile .Rprofile

# Install renv and restore the project library ---
# RUN Rscript -e 'renv::restore(lockfile = "renv.lock", repos = NULL)'
# RUN Rscript -e 'install.packages(c("renv", "BiocManager"))'

# Install Source Sans 3 (https://fonts.google.com/specimen/Source+Sans+3)
# RUN mkdir -p /usr/share/fonts/truetype/google-source-sans-3 && \
#     cd /usr/share/fonts/truetype/google-source-sans-3 && \
#     curl -L -O https://github.com/adobe-fonts/source-sans/archive/refs/tags/3.052R.tar.gz && \
#     tar -xzf 3.052R.tar.gz --strip-components=1 && \
#     rm 3.052R.tar.gz
# RUN fc-cache -f -v
RUN R -e "install.packages(c('extrafont', 'sysfonts', 'showtext', 'remotes'))"
#RUN R -e "extrafont::font_import(prompt = FALSE, pattern = 'SourceSans')"

# Copy the rest of the pipeline code ---
COPY utils/ utils/
COPY gwas-make-plots.R gwas-make-plots.R
COPY effect-size-comparison.R effect-size-comparison.R