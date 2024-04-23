### In this file we generate the datasets for the linear Pareto test ###

# Namely we use the centroids of the clusters from scMultiSim as the archetypes
# And generate the by sampling from a simplex with the archetypes as vertices

library(scMultiSim)
library(tidyverse)
library(foreach)
library(doParallel)
data("GRN_params_100")
source("utils.R")

registerDoParallel(3)


data_dir <- "../data/phyla3_no_batch_pareto"

dir.create(data_dir, showWarnings = F)

N_sims <- 20

noises <- c(0.6,0.75,0.9)

# the actual simulation
options <- list(
  rand.seed = 3,
  GRN = GRN_params_100,
  num.cells = 1000,
  num.cifs = 100,
  cif.sigma = 0.5,
  tree = Phyla5(),
  diff.cif.fraction = 0.8,
  do.velocity = F,
  discrete.cif = T
  
)


# phyla 3 simulations

for(i in 1:N_sims) {
  x <- foreach(j = noises) %dopar% {
    set.seed(i)
    nms <- paste0("dataset", i)
    tree <- Phyla3()
    results <- sim_true_counts(
      options %>% purrr::list_modify(diff.cif.fraction = j, tree = tree, seed = i, rand.seed = i))
    results$tree <- tree
    
    # We get RNA and ATAC centroids
    
    dimnames(results$atac_counts) <- list(paste0("peak", 1:nrow(results$atac_counts)), results$cell_meta$cell_id)
    dimnames(results$counts) <- list(paste0("gene", 1:nrow(results$counts)), results$cell_meta$cell_id)
    archetypes_ATAC <- sapply(unique(results$cell_meta$pop), groupwise_means, df=results$cell_meta, mat=results$atac_counts)
    colnames(archetypes_ATAC) <- paste0("arch", 1:ncol(archetypes_ATAC))
    archetypes_RNA <- sapply(unique(results$cell_meta$pop), groupwise_means, df=results$cell_meta, mat=results$counts, round = TRUE)
    colnames(archetypes_RNA) <- paste0("arch", 1:ncol(archetypes_RNA))
    
    narchetypes <- ncol(archetypes_RNA)
    archetypes <- list(RNA = archetypes_RNA, ATAC = archetypes_ATAC) 
    
    # Here we sample from a Dirichlet the weights
    
    results$archetypes <- archetypes
    A <- gtools::rdirichlet(options$num.cells, rep(1, narchetypes))
    results$A_true <- A
    results$cell_meta$pop <- paste0("arch", apply(A,1, which.max))
    
    # We reconstruct the simplex by doing A %*% archetypes
    
    results$counts <- archetypes_RNA %*% t(A)
    results$atac_counts <- archetypes_ATAC %*% t(A)
    
    expr_dnames <- dimnames(results$counts)
    results$counts <- sapply( 1:ncol(results$counts), function(i) rpois(length(results$counts[,i]), results$counts[,i]) )
    dimnames(results$counts)  <- expr_dnames
    
    atac_dnames <- dimnames(results$atac_counts)
    results$atac_counts <- sapply( 1:ncol(results$atac_counts), function(i) rpois(length(results$atac_counts[,i]), results$atac_counts[,i]) )
    dimnames(results$atac_counts)  <- atac_dnames
    
    # We add expressoion noise
    
    add_expr_noise(results)

    dir_tmp <- file.path(data_dir, nms)
    dir.create(dir_tmp, showWarnings = FALSE, recursive = T)
    file_name <- file.path(dir_tmp, paste0("simulation_result_cif_",j, ".rds"))
    saveRDS(results, file = file_name) 
  }
}


