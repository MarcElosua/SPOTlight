FROM rocker/rstudio:3.6.3

# RUN mkdir /home/analysis

#### Defining arguments to pass ####
# ARG WHEN

#### Install SpaceDecon packages ####
# libcurl4-openssl-dev
RUN apt-get update -y && \
 apt-get install -y libcurl4-gnutls-dev \
libssl-dev \
libssh2-1-dev \
libxml2-dev \
libhdf5-dev \
libgmp-dev \
libpng-dev \
libgsl-dev \
libxt-dev \
libcairo2-dev \
libtiff-dev \
fftw-dev

#### Install CRAN packages ####
RUN R -e "options(repos = \
 list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2020-09-15')); \
 install.packages(c('arrangements', 'cowplot', 'doSNOW', 'doParallel', \
 'dplyr', 'dtplyr', 'foreach', 'ggplot2', 'ggpubr', 'igraph', 'imager', \
 'matrixStats', 'knitrProgressBar', 'Matrix', 'NMF', 'nnls', 'pdist', \
 'philentropy', 'purrr', 'Seurat', 'tibble', 'httr', \
 'metap', 'plotly', 'git2r', 'png', 'hdf5r'))"

#### Install bioconductor packages ####
RUN R -e "options(repos = \
  list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2020-09-15')); \
  install.packages('BiocManager')"

#### Install bioconductor packages ####
RUN R -e "options(repos = \
  list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2020-09-15')); \
  BiocManager::install(pkgs = c('SingleCellExperiment', 'BiocGenerics', \
 'DropletUtils', 'multtest', 'edgeR'), version = '3.10')"

RUN R -e "options(repos = \
 list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2020-09-15'), \
      BiocManager::repositories()); \
 install.packages(c('git2r', 'usethis' , 'devtools', 'testthat', 'roxygen2')); \
 devtools::install_github('https://github.com/MarcElosua/SPOTlight')"
