library(scMultiSim)
library(tidyverse)
library(foreach)
library(doParallel)
data("GRN_params_100")


registerDoParallel(3)

N_sims <- 10

n_CELLS <- c(1000, 5000, 10000, 25000, 50000, 100000)

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


my_phyla_4 <- function(){
  library(dplyr)
  tree <- list()
  tree$edge <- data.frame(c(5,5,5,5), c(1,2,3,4)) %>% as.matrix()
  tree$Nnode <- 1
  tree$edge.length <- c(1,1,1,1)
  tree$tip.label <- c(1,2,3,4)
  class(tree) <- "phylo"
  return(tree)
}


data_dir <- "../data_timing/"

dir.create(data_dir, showWarnings = F)

# timing simulations

for(i in 1:N_sims) {
  for(j in n_CELLS){
    set.seed(i)
    nms <- paste0("dataset", i)
    tree <- my_phyla_4()
    results <- sim_true_counts(
      options %>% purrr::list_modify(diff.cif.fraction = 0.75, tree = tree, seed = i, rand.seed = i))
    results$tree <- tree
    
    # expr_dnames <- dimnames(results$counts)
    # results$counts <- sapply( 1:ncol(results$counts), function(i) rpois(length(results$counts[,i]), results$counts[,i]) )
    # dimnames(results$counts)  <- expr_dnames
    # 
    # atac_dnames <- dimnames(results$atac_counts)
    # results$atac_counts <- sapply( 1:ncol(results$atac_counts), function(i) rpois(length(results$atac_counts[,i]), results$atac_counts[,i]) )
    # dimnames(results$atac_counts)  <- atac_dnames
    
    results$counts <- log1p(results$counts)
    results$atac_counts <- log1p(results$atac_counts)
    
    dir_tmp <- file.path(data_dir, nms)
    dir.create(dir_tmp, showWarnings = FALSE, recursive = T)
    file_name <- file.path(dir_tmp, paste0("simulation_result_ncells_",j,"_rep_",i ,".rds"))
    saveRDS(results, file = file_name) 
  }
}
