# Transdimensional Predictive Least Squares (TPLS)

Note: Python implementation at [tpls_python](https://github.com/marija-iloska/tpls_python).

This code is a MATLAB implementation of the algorithm TPLS proposed in our paper "Transdimensional Model Learning with Online Feature Selection based on Predictive Least Squares".
We provide an example code on how a user can run TPLS, as well as scripts to reproduce the results and figures in the paper.

## Introduction
TPLS is a distribution-free online feature selection algorithm that is completely based on least squares (LS). With new data arrivals, TPLS recursively updates the parameter estimate not only
in value, but in dimension as well. What makes TPLS unique is the ability to recursively move up and down model dimension, and the fact that it uses the predictive (instead of the fitting) error as its criterium whether to add or remove features.
Specifically, the foundations of TPLS are recursive LS [(RLS)](https://dl.acm.org/doi/book/10.5555/151045), order recursive LS [(ORLS)](https://dl.acm.org/doi/book/10.5555/151045) and predictive LS [(PLS)](https://academic.oup.com/imamci/article-abstract/3/2-3/211/660741).

## Code

### How to run TPLS
Required files: <br/> 
example_code.m <br/> 
util/ <br/> 
tpls/ <br/> 
<br/> 
Script to run: <br/>
example_code.m - a script that demonstrates how to call and run TPLS. It includes feature bar plots and predictive error. <br/> 
<br/> 

### jpls/
ls_updates/ - a folder that contains the least squares update algorithms: <br/> 
time_update.m  - implements RLS (update with arrival of new data point) <br/> 
ols_ascend.m  - implements ascending ORLS (increase model dimension by 1 --> add 1 feature) <br/> 
ols_descend.m - implements descending ORLS (decrease model dimension by 1 --> remove 1 feature) <br/> 

jumps_up.m  and  jumps_down.m - functions which find the predictive error for adding/removing each feature to/from the current model. <br/> 
pred_error.m - a function that computes the predictive error recursively in dimension. <br/> 
expectations.m - a function that computes the MSE at single time instants for [regret analysis](https://pubsonline.informs.org/doi/abs/10.1287/opre.30.5.961). <br/> 

### Reproduce Results
Scripts to run: <br/>
reproduce_fig1_fig4_fig5.m <br/>
reproduce_fig2_fig3.m <br/>

These scripts run TPLS against the baseline methods [OLinLASSO](https://proceedings.mlr.press/v206/yang23g.html) and [rjMCMC](https://ieeexplore.ieee.org/abstract/document/7089644?casa_token=SiWa6_nyqegAAAAA:Kmx-raPBJg8OPWSNuAT5UXbUtAQ5DXdzmgbg2N8lm2JCkKIlIvLNDMY4AE_Bc80_FU8wAFylng), as well as the ground truth. The results in figs 1,4,5 are single runs, where as those in figs 2 and 3 are statistical, averaged over R=100 runs. <br\>




