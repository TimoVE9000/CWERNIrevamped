
# CWERNIrevamped

<!-- badges: start -->
[![R-CMD-check](https://github.com/TimoVE9000/CWERNIrevamped/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TimoVE9000/CWERNIrevamped/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
The original [CWERNI package](https://github.com/TimoVE9000/CWERNI) was developed to simulate community wide evolutionary rescue in a neutral community. It accompanied the publication: Van Eldijk, Bisschop & Etienne (2020), Uniting Community Ecology and Evolutionary Rescue Theory: Community-Wide Rescue Leads to a Rapid Loss of Rare Species, Frontiers in Ecology and Evolution, Volume 8 - 2020 [DOI](https://doi.org/10.3389/fevo.2020.552268). The goal of CWERNIrevamped is to re-make the the original package, making it compatible with R 4.5 and using best development pratices. This includes writing good documentation, proper use of git & github and adding some automated testing. 


## Installation

You can install the development version of CWERNIrevamped from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("TimoVE9000/CWERNIrevamped")
```

## Example

This example shows how to simulate community wide rescue in a neutral community. For documentation of the other functions in the package see the [manual](https://github.com/TimoVE9000/CWERNIrevamped/manual.pdf)

``` r
library(CWERNIrevamped)

##Generate a neutral community
Community = generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000))

##Simulate 



```

