
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SUFA

A `C++`-based implementation of the SUbspace Factor Analysis (SUFA)
model by Chandra, Dunson, and Xu (2023).

<!-- badges: start -->
<!-- badges: end -->

## Installation

The package has been tried and tested in Ubuntu and macOS.

### Required System Packages

- `openmp`

- `R (>= 4.3.1)`

### Installation from Github

You can install the development version of SUFA from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("noirritchandra/SUFA", build_vignettes = TRUE)
```

Must set `build_vignettes = TRUE` to install the vignette files
containing illustrations. However, this may take considerably more time
to install the package.

## Simulation Example

Refer to the following vignette for illustrations in simulated examples:

``` r
vignette(topic="Intro_simulation",package = "SUFA")
```

## Genedata Application

Refer to the following vignette in gene expressions datasets for
inference on gene networks:

``` r
vignette(topic="Genedata_application",package = "SUFA")
```

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-chandra2023sufa" class="csl-entry">

Chandra, Noirrit Kiran, David B Dunson, and Jason Xu. 2023. “Inferring
Covariance Structure from Multiple Data Sources via Subspace Factor
Analysis.” *arXiv:2305.04113*.

</div>

</div>
