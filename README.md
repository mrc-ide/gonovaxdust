# gonovaxdust

<!-- badges: start -->
[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![R build status](https://github.com/mrc-ide/gonovaxdust/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/gonovaxdust/actions)
<!-- badges: end -->

This package implements a stochastic compartmental model of gonorrhoea infection with vaccination.

## Installation

You will need a compiler to install dependencies for the package, and to build
the models. Use `pkgbuild::check_build_tools()` to see if your system is usable.

You will need the packages `odin.dust`, which can be installed using:

```r
remotes::install_github("mrc-ide/odin.dust", upgrade = FALSE)
```


The package can then be installed directly from GitHub with:

```r
remotes::install_github("mrc-ide/gonovaxdust", upgrade = FALSE)
```

## License

MIT © Imperial College of Science, Technology and Medicine

