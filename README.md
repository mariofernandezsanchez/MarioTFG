# MarioTFG

This package is made by Mario Fernández Sánchez for the achievement of his TFG. Whose purpose is to standardize the process of prediction of the most likely motif, develop a scoring system for each peaks, and determine a threshold to separate the peaks potentially functional of those that are not. From work previously done by: Ginés Almagro-Hernández, Jesualdo Tomás Fernández-Breis EMAILS: [gines.almagro\@um.es](mailto:gines.almagro@um.es){.email}, [jfernand\@um.es](mailto:jfernand@um.es){.email}

## Installation

En primer lugar, tendremos que instalar el paquete `devtools` como hemos visto en un apartado anterior:

``` r
install.packages("devtools")
```

Once we have `devtools` installed, we will use the functions for installing other packages. The options are as follows:

- `install_bioc()` from [Bioconductor](https://www.bioconductor.org/),
- `install_bitbucket()` from [Bitbucket](https://bitbucket.org/),
- `install_cran()` from [CRAN](https://cran.r-project.org/index.html),
- `install_git()`from a [git](https://git-scm.com/) repository,
- `install_github()` from [GitHub](https://github.com/),
- `install_local()` from a file hosted on our computer,
- `install_svn()` from an [SVN](https://subversion.apache.org/) repository,
- `install_url()` from a URL, and
- `install_version()` for a specific version of a CRAN package.

For example, to install the [MarioTFG package](https://github.com/mariofernandezsanchez/MarioTFG) from its github repository we can do it as shown below:

``` r
devtools::install_github("mariofernandezsanchez/MarioTFG")
```

## R

In this directory is the file functions_r, which is where all the functions that are part of the package are found.

## EtapasTFG

In this directory In this directory you will find the results obtained in the realization of the TFG
