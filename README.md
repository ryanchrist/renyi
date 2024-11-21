# renyi
This package implements the Renyi Outlier Test [arXiv:2411.13542](https://arxiv.org/abs/2411.13542), designed for modern large scale testing applications, especially those where prior information available.

The test combines a vector of independent uniform p-values into one p-value with power 
against alternatives where a small number of p-values are non-null.

The test can leverage prior probabilities/weights specifying which variables are likely to be outliers and prior estimates of effect size. 
The procedure is fast even when the number of initial p-values is large (e.g. in the millions) and numerically stable even for very small p-values (e.g. 10^-300).
