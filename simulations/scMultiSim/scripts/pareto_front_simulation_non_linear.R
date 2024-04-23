### In this file we generate the datasets for the non-linear Pareto test ###

# We use the centroids of the clusters from scMultiSim after projecting them to a non-linear latent space using a VAE 
# as the archetypes and generate the by sampling from a simplex with the archetypes as vertices

library(scMultiSim)
library(tidyverse)
library(foreach)
library(doParallel)
library(kernlab)
data("GRN_params_100")
source("utils.R")

# You can use the environment defined in the requirements.txt here

reticulate::use_condaenv("scdeepaa_sim", required = T)
sys <- reticulate::import("sys")
sys$path <- c(sys$path, "/home/salvatore.milite/work/python_packages/multiDeepAA/src/")

scdeepaa <- reticulate::import("scdeepaa")
torch <- reticulate::import("torch")


registerDoParallel(1)


data_dir <- "../data/phyla3_no_batch_nlpareto"

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
    
    dimnames(results$atac_counts) <- list(paste0("Peak", 1:nrow(results$atac_counts)), results$cell_meta$cell_id)
    dimnames(results$counts) <- list(paste0("Gene", 1:nrow(results$counts)), results$cell_meta$cell_id)
    
    # We scale teh data as it works better with VAEs but we keep
    # means and variances of the features to reconstruct the normalized counts back
    
    views <- list(rna = results$counts  %>% t %>% scale, 
                  atac = results$atac_counts %>% t %>% scale)
    
    rownames(views$rna) <- paste0("Cell", 1:nrow(views$rna)) 
    rownames(views$atac) <- paste0("Cell", 1:nrow(views$atac)) 
    
    colnames(views$rna) <- paste0("Gene", 1:ncol(views$rna)) 
    colnames(views$atac) <- paste0("Peak", 1:ncol(views$atac)) 
    
    n_archetypes <- length(unique(results$cell_meta$pop))
    
    # We fit a VAE as a non-linear function
    
    MULT = 16L
    
    res_auto_hidden_dims <- scdeepaa$fit_deepAA(
      list(views$rna, views$atac),
      list(rep(1., dim(views$rna)[1]), rep(1., dim(views$atac)[1])),
      list("G", "G"),
      hidden_dims_dec_common = list(8L * MULT,16L * MULT),
      hidden_dims_dec_last = list(32L * MULT ,64L * MULT),
      hidden_dims_enc_ind = list(64L * MULT, 32L * MULT),
      hidden_dims_enc_common = list(16L * MULT),
      hidden_dims_enc_pre_Z = list(8L * MULT),
      lr = 0.001,
      gamma_lr = 0.1,
      steps = 1000L,
      narchetypes = as.integer(n_archetypes),
      batch_size = 9000L,
      fix_Z = FALSE,
      just_VAE = TRUE
    )
    
    Z_all <- res_auto_hidden_dims$inferred_quantities$Z %>%  t
    colnames(Z_all) <-paste0("cell", 1:nrow(views$rna)) 
    
    archetypes_Z <- sapply(unique(results$cell_meta$pop), groupwise_means, df=results$cell_meta, 
                           mat=Z_all)
    
    archetypes <- archetypes_Z
    
    # Generate the archetypes weights
    
    results$archetypes <- archetypes
    A <- gtools::rdirichlet(options$num.cells, rep(1, n_archetypes))
    results$A_true <- A
    results$cell_meta$pop <- paste0("arch", apply(A,1, which.max))
    
    # We reconstruct the counts from the VAEs fit
    
    reconstructed_arcs <- res_auto_hidden_dims$deepAA_obj$decoder(torch$tensor(t(archetypes_Z %*% t(A)))$float())
    
    archetypes_RNA <- reconstructed_arcs[[1]][[1]]$detach()$cpu()$numpy()
    archetypes_ATAC <- reconstructed_arcs[[1]][[2]]$detach()$cpu()$numpy()
    
    # We go from scaled data to normalized counts again
    
    archetypes_RNA <- sapply(1:ncol(archetypes_RNA), function(i) {
      (archetypes_RNA[,i] * attributes(views$rna)[[4]][i]) + attributes(views$rna)[[3]][i]
    })
    
    archetypes_ATAC <- sapply(1:ncol(archetypes_ATAC), function(i) {
      (archetypes_ATAC[,i] * attributes(views$atac)[[4]][i]) + attributes(views$atac)[[3]][i]
    })
    
    rownames(archetypes_ATAC) <- colnames(results$atac_counts)
    
    colnames(archetypes_ATAC) <- paste0("arch", 1:ncol(archetypes_ATAC))
    
    rownames(archetypes_RNA) <- colnames(results$counts)
    
    colnames(archetypes_RNA) <- paste0("arch", 1:ncol(archetypes_RNA))
    
    
    narchetypes <- ncol(archetypes_RNA)

    results$counts <- archetypes_RNA %>% t
    results$atac_counts <- archetypes_ATAC %>% t
    
    # We add some Poisson sequencing noise
    
    expr_dnames <- dimnames(results$counts)
    results$counts <- sapply( 1:ncol(results$counts), function(i) rpois(length(results$counts[,i]), results$counts[,i]) )
    dimnames(results$counts)  <- expr_dnames
    
    atac_dnames <- dimnames(results$atac_counts)
    results$atac_counts <- sapply( 1:ncol(results$atac_counts), function(i) rpois(length(results$atac_counts[,i]), results$atac_counts[,i]) )
    dimnames(results$atac_counts)  <- atac_dnames
    
    results$counts_obs <- results$counts
    results$atacseq_obs <- results$atac_counts
    
    dir_tmp <- file.path(data_dir, nms)
    dir.create(dir_tmp, showWarnings = FALSE, recursive = T)
    file_name <- file.path(dir_tmp, paste0("simulation_result_cif_",j, ".rds"))
    saveRDS(results, file = file_name) 
  }
}


