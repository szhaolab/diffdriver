# diffDriver

## To install

``` remotes::install_github("szhaolab/diffdriver",ref ="main") 

```

## How to use 

Please visit the following website for the
details of this pacakge, [diffdriver manual and
vignette](https://szhaolab.github.io/diffdriver/).  We have
included example scripts to run `diffdriver`.  Please take a
look at `scripts/` folder. 

* `run_nttype.R` gives an example of how to run `diffdriver`
  under regular model, where the annotation files are
categorized into nine catogories based on its DNA mutation type. Another script  `run_signature.R` givens an example
of how to run `diffdriver` under sinature mode.You can run the scripts as
`Rscript run_diffdriver.R.  Note, we
suggested to run drivermaps in your actual run in order to
prepare the required inputs for diffdriver.

* `simulate_functions.R` gives the functions for simulation studies of `diffdriver`.  We generate mutation and phenotype
data for a single gene many times, assuming it is under
differential selection or not, and check how many times we
can detect it using. Three scenarioes are considered in the simulations, one for phenotypes which are independent of backgound mutation, 
two for phenotypes which are correlated with background mutation. 

