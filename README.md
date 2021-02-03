<img src="img/SPOTlight_VF2.png" width="200px" style="display: block; margin: auto;" />

The goal of **SPOTlight** is to provide a tool that enables the
deconvolution of cell types and cell type proportions present within
each capture locations comprising mixtures of cells, originally
developed for 10X’s Visium - spatial trancsiptomics- technology, it can
be used for all technologies returning mixtures of cells. **SPOTlight**
is based on finding topic profile signatures, by means of an NMFreg
model, for each cell type and finding which combination fits best the
spot we want to deconvolute.

![Graphical abstract](img/SPOTlight_graphical_abstract.png)

Installation
------------

You can install the latest stable version from the GitHub repository
[SPOTlight](https://github.com/MarcElosua/SPOTlight) with:

``` r
# install.packages("devtools")
devtools::install_github("https://github.com/MarcElosua/SPOTlight")
devtools::install_git("https://github.com/MarcElosua/SPOTlight")
```

Or the latest version in development by downloading the devel branch

``` r
devtools::install_github("https://github.com/MarcElosua/SPOTlight", ref = "devel")
devtools::install_git("https://github.com/MarcElosua/SPOTlight", ref = "devel")
```

To ensure the environment is compatible we have put out a docker
environment that can be downloaded from DockerHub.  
To download the R environment image

    # R environment
    docker pull marcelosua/spotlight_env_r:latest

To download the R studio image

    # Rstudio environment
    docker pull marcelosua/spotlight_env_rstudio:latest

Run the following command in the terminal

    docker run -e PASSWORD=pwd -p 8787:8787 marcelosua/spotlight_env_rstudio 

Go to
**<a href="http://localhost:8787/" class="uri">http://localhost:8787/</a>**
or other port you’ve set before.

Log in with:

**username**: rstudio

**password**: pwd \# Same password as set in the above command (pwd)

In the Rstudio environment you can upload files and carry out analysis
there, to save the files from the RStudio server click on More (gear)
and export the desired files.

For more information check out the
[rocker](https://www.rocker-project.org/) project guidelines
<a href="https://hub.docker.com/r/rocker/rstudio/" class="uri">https://hub.docker.com/r/rocker/rstudio/</a>

A vignette on how to run the basic SPOTlight workflow can be found [here](https://marcelosua.github.io/)
