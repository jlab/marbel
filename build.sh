#!/bin/bash

$PYTHON -m pip install . 

Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); BiocManager::install('polyester')"
