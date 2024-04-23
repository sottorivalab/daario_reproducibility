# Synthetic simulation 

In this folder we have the scripts to reproduce the analysis on simulated data. This part is in R as our simulation tool `scMultiSim` is R based, we will call MIDAA using reticulate. As such to run successfully this piece you need the general environment defined in the root directory of this repository but also an R session with the following packages installed:

```r
install.packages("reticulate", "ddevtools", "tidyverse", 
      "RColorBrewer", "matrixStats", "foreach",
      "doParallel", "scales", "kernlab")
      
devtools::install_github("ZhangLabGT/scMultiSim")
```

After you have installed the packages, you first have to run the 3 scripts to generate the synthetic counts:
* `pareto_front_simulations_linear.R`: generates 
* `pareto_front_simulations_non_linear.R`: generates 
* `branching_simulations.R`: generates 

Once you have done this 
