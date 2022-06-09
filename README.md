# diffDriver

## To install

```
remotes::install_github("szhaolab/diffdriver",ref = "main")
```


## How to use

We have included example scripts to run `diffdriver`. Please take a look at `scripts/` folder. 

* `run_diffdriver.R` gives an example of how to run `diffdriver` using example files provided with the package. You can download this script and run it as `Rscript run_diffdriver.R` in command line. Note, we suggested to run drivermaps in your actual run in order to prepare the required inputs for diffdriver.

* `run_simulate_singlegene.R` gives an example of how to perform simulations. We generate mutation and phenotype data for a single gene many times, assuming it is under differential selection or not, and check how many times we can detect it using. You can download this script and run it as `Rscript run_simulate_singlegene.R` in command line. 

