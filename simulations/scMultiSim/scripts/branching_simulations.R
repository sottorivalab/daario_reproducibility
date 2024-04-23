### In this file we generate the datasets for the Trajectory inference case ###

# Here luckily scMultiSim does the simulation for us
# We just specify the tree structure, namely 3 branches of differentiation 
# From an unique stem population

library(scMultiSim)
library(tidyverse)
library(foreach)
library(doParallel)
data("GRN_params_100")


registerDoParallel(3)


data_dir <- "../data/phyla3_no_batch_trajectory"

dir.create(data_dir, showWarnings = F)

N_sims <- 20

noises <- c(0.6,0.75,0.9)

options <- list(
  rand.seed = 3,
  GRN = GRN_params_100,
  num.cells = 1000,
  num.cifs = 100,
  cif.sigma = 0.5,
  tree = Phyla5(),
  diff.cif.fraction = 0.8,
  do.velocity = T
)




my_phyla_3 <- function(){
  library(dplyr)
  tree <- list()
  tree$edge <- data.frame(c(5,5,5), c(1,2,3)) %>% as.matrix()
  tree$Nnode <- 1
  tree$edge.length <- c(1,1,1)
  tree$tip.label <- c(1,2,3)
  class(tree) <- "phylo"
  return(tree)
}



# phyla 3 simulations

for(i in 1:N_sims) {
  x <- foreach(j = noises) %dopar% {
    set.seed(i)
    nms <- paste0("dataset", i)
    tree <- my_phyla_3()
    results <- sim_true_counts(
      options %>% purrr::list_modify(diff.cif.fraction = j, tree = tree, seed = i, rand.seed = i))
    results$tree <- tree
    
    results$counts <- log1p(results$counts)
    results$atac_counts <- log1p(results$atac_counts)
    
    dir_tmp <- file.path(data_dir, nms)
    dir.create(dir_tmp, showWarnings = FALSE, recursive = T)
    file_name <- file.path(dir_tmp, paste0("simulation_result_cif_",j, ".rds"))
    saveRDS(results, file = file_name) 
  }
}



