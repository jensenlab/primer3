
primer3
=======

primer3 provides a direct interface to the thermodynamic calculations in Primer3 (melting temperatures, primer dimers, secondary structures, etc.). The package is self-contained and does not require a separate installation of Primer3.

Overview
--------

primer3 contains four functions:

``` r
calculate_tm()
calculate_hairpin()
calculate_homodimer()
calculate_dimer()
```

See the function documentation for details on inputs and usage.

Installation
------------

``` r
# Install the latest version directly from GitHub
install.packages("devtools")
devtools::install_github("jensenlab/primer3")
```

Usage
-----

``` r
calculate_tm(c("AAGCCGCGTACGA", "AAGAGCGATGACG"))
## [1] 44.8188 38.7431

calculate_hairpin("CCCCCATCCGATCAGGGGG")
## $structure_found
## [1] TRUE
## 
## $temp
## [1] 62.62755
## 
## $ds
## [1] -108.1073
## 
## $dh
## [1] -36300
## 
## $dg
## [1] -2770.524
## 
## $align_end_1
## [1] -36300
## 
## $align_end_2
## [1] -104
```

Acknowledgements
----------------

primer3 builds on [Primer3](http://primer3.sourceforge.net) and modifications made by the [primer3-py](https://github.com/libnano/primer3-py) project.
