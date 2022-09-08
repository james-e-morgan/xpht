
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Extended Persistent Homology Transform

<!-- badges: start -->
<!-- badges: end -->

XPHT allows the computation of the extended persistent homology
transform of binary images detailed in [Turner, Robins, & Morgan
(2022)](doi:10.48550/arXiv.2208.14583). This provides a simple method of
computing the Wasserstein $p$ distance between image for statistical
shape analysis.

## Installation

This package relies on the [Imager
package](https://cran.r-project.org/package=imager). This allows us to
read in images in order to extract the boundary curves. There is
detailed information about installing this package, so we direct the
user to view that if any difficulties are encountered.

If you wish to run the development version, install the devtools package
(if you haven’t already) and then run:

``` r
devtools::install_github("james-e-morgan/xpht")
```

## A Note on Reading Images

XPHT requires the input be binary images. For this it is required that
the pixel values of the input images are either 0 or 1. If any other
values are present an error will be thrown.
