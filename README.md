
The dccpp R Package
======================

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/github/actions/workflow/status/berrij/dccpp/R-CMD-check.yaml?&style=for-the-badge)](https://github.com/BerriJ/dccpp/actions/workflows/R-CMD-check.yaml)
[![GitHub Workflow Status (branch)](https://img.shields.io/github/actions/workflow/status/berrij/dccpp/pkgdown.yaml?label=Documentation&style=for-the-badge)](https://dccpp.berrisch.biz/)
[![Lifecycle: stable](https://img.shields.io/badge/Lifecycle-stable-brightgreen?style=for-the-badge)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

Fast computation of the distance covariance `dcov` and distance correlation `dcor`. The computation cost is only O(n log(n)) for the distance correlation (see Chaudhuri, Hu, 2019, [arXiv](https://arxiv.org/abs/1810.11332), [elsevier](https://doi.org/10.1016/j.csda.2019.01.016)). The functions are written entirely in C++ to speed up the computation.


Installation
------------

### Install from CRAN

``` r
install.packages("dccpp")
```

### Install from GitHub

``` r
# install.packages("remotes")
remotes::install_github("BerriJ/dccpp")
```

## Usage

``` r
dcov(x, y) # To calculate the distance covariance
dcor(x, y) # To calculate the distance correlation
```

## Contributions and Issues

Feel free to [raise an issue](https://github.com/BerriJ/dccpp/issues/new) if you find something not working properly.

You are also very welcome to contribute to `dcccp`. Please base your pull requests on the development branch. Note that this package focuses on performance, PR's that improve the performance are particularly welcome.

## License

[GNU General Public License](https://www.gnu.org/licenses/) (≥ 3)
