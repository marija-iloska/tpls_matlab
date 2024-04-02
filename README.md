# Jump Predictive Least Squares (JPLS)

Note: Python implementation at [JPLS_Python](https://github.com/marija-iloska/JPLS).

This code is a MATLAB implementation of the algorithm JPLS proposed in our paper "Transdimensional Model Learning with Online Feature Selection based on Predictive Least Squares".
We provide an example code on how a user can run JPLS, as well as scripts to reproduce the results and figures in the paper.

## Introduction
JPLS is a distribution-free online feature selection algorithm that is completely based on least squares (LS). With new data arrivals, JPLS recursively updates the parameter estimate not only
in value, but in dimension as well. What makes JPLS unique is the ability to recursively jump up and down model dimension, and the fact that it uses the predictive (instead of the fitting) error as its criterium whether to add or remove features.
Specifically, the foundations of JPLS are recursive LS [(RLS)](https://dl.acm.org/doi/book/10.5555/151045), order recursive LS [(ORLS)](https://dl.acm.org/doi/book/10.5555/151045) and predictive LS [(PLS)](https://academic.oup.com/imamci/article-abstract/3/2-3/211/660741).

## Code

