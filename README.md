# Jump Predictive Least Squares (JPLS)

Note: Python implementation at [jpls_python](https://github.com/marija-iloska/JPLS).

This code is a MATLAB implementation of the algorithm JPLS proposed in our paper "Transdimensional Model Learning with Online Feature Selection based on Predictive Least Squares".
We provide an example code on how a user can run JPLS, as well as scripts to reproduce the results and figures in the paper.

## Introduction
JPLS is a distribution-free online feature selection algorithm that is completely based on least squares (LS). With new data arrivals, JPLS recursively updates the parameter estimate not only
in value, but in dimension as well. What makes JPLS unique is the ability to recursively jump up and down model dimension, and the fact that it uses the predictive (instead of the fitting) error as its criterium whether to add or remove features.
Specifically, the foundations of JPLS are recursive LS [(RLS)](https://dl.acm.org/doi/book/10.5555/151045), order recursive LS [(ORLS)](https://dl.acm.org/doi/book/10.5555/151045) and predictive LS [(PLS)](https://academic.oup.com/imamci/article-abstract/3/2-3/211/660741).

## Code

### How to run JPLS
Required files: <br/> 
example_code.m <br/> 
util/ <br/> 
jpls/ <br/> 
<br/> 
example_code.m - a script that demonstrates how to call and run JPLS. It includes feature bar plots and predictive error. <br/> 
<br/> 

### jpls/
ls_updates/ - a folder that contains the least squares update algorithms: <br/> 
time_update.m  - implements RLS (update with arrival of new data point) <br/> 
orls_ascend.m  - implements ascending ORLS (increase model dimension by 1 --> add 1 feature) <br/> 
orls_descend.m - implements descending ORLS (decrease model dimension by 1 --> remove 1 feature) <br/> 
<br/> 

jumps_up.m  and  jumps_down.m - functions which find the predictive error for adding/removing each feature to/from the current model. <br/> 
pred_error.m - a function that computes the predictive error. <br/> 
expectations.m - a function that computes the MSE at single time instants for [regret analysis](https://pubsonline.informs.org/doi/abs/10.1287/opre.30.5.961). <br/> 

###


