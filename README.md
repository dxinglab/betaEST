# betaEST
This package provide R functions for the following purposes: 

`betax`: for the estimation of the new unbiased multi-community beta index. the input (named `comp`) must be a matrix or data frame of species abundances with species on the rows and sites on the columns.

`betap`: estimating Whittaker's classic indices (both gamma/alpha, returned with name betaw, and 1-alpha/gamma, returned with name betap). the input (named `comp`) must be a matrix or data frame of either species abundances or incidences with species on the rows and sites on the columns.

`betadev`: estimating betadev from 0/1 data. there are three inputs: `comp` must be a matrix or data frame of species incidences with species on the rows and sites on the columns. this can also be species abundances, but the abundance matrix will be reduced to incidence matrix. `S`: a single number for total species richness in the metacommunity. `N`: a single number for total species abundance in the metacommunity.
