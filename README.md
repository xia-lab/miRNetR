# miRNetR: A companion R package for the miRNet web server

Note - miRNetR is still under development - we cannot guarantee full functionality 

## Description

**miRNetR** contains the R functions and libraries underlying the popular miRNet web server. After installing and loading the package, users will be able to reproduce the same results from their local computers using the corresponding R command history downloaded from miRNet, thereby achieving maximum flexibility and reproducibility.

## Getting Started

### Step 1. Install package dependencies

To use miRNetR , first install all package dependencies. Ensure that you are able to download packages from bioconductor. To install package dependencies, use the pacman R package (for those with >R 3.5.1). Note that some of these packages may require additional library dependencies that need to be installed prior to their own successful installation.

```R
install.packages("pacman")

library(pacman)

pacman::p_load(RSQLite, Cairo, fastmatch, igraph, RJSONIO, foreach, doParallel, preprocessCore, limma, edgeR, HTqPCR, genefilter)
```
### Step 2. Install the package

miRNetR is freely available from GitHub. The package documentation, including the vignettes for each module and user manual is available within the downloaded R package file. If all package dependencies were installed, you will be able to install the miRNetR. Due to issues with Latex, some users may find that they are only able to install miRNetR without any documentation (i.e. vignettes).

Install the package directly from github using the *devtools* package. Open R and enter:

```R
# Step 1: Install devtools
install.packages("devtools")
library(devtools)

# Step 2: Install miRNetR WITHOUT documentation
devtools::install_github("xia-lab/miRNetR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"))

# Step 2: Install miRNetR WITH documentation
devtools::install_github("xia-lab/miRNetR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)

```
## Tutorials

For detailed tutorials on how to use miRNetR, please refer to the R package vignettes. Note, the functions below work only if the R package vignettes were built.

Within R:
```R
vignette(package="miRNetR")
```

Within a web-browser:
```R
browseVignettes("miRNetR")
```

## Citation

miRNetR has been developed by the [XiaLab](http://xialab.ca/) at McGill University. The original manuscript (web-based version) can be found [here](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkaa467/5850315).

We encourage users to further develop the package to suit their needs. If you use the R package, please cite us:

Chang, L., Zhou, G., Soufan, O. and Xia, J. (2020) miRNet 2.0 - network-based visual analytics for miRNA functional analysis and systems biology. Nucl. Acids Res. (doi: 10.1093/nar/gkaa467)

Fan Y, Siklenka, K., Arora, SK., Ribeiro, P., Kimmins, S. and Xia, J. (2016) miRNet - dissecting miRNA-target interactions and functional associations through network-based visual analysis. Nucl. Acids Res. 44 W135–141

## Bugs or feature requests

To inform us of any bugs or requests, please open a new issue or send an email to #le.chang@mail.mcgill.ca.

## To Do

mRNA-miRNA co-expression network analysis 

Network motif analysis 
