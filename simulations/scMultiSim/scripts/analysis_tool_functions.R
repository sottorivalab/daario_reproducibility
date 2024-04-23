library(dplyr)

# dir = "/home/salvatore.milite/data/multideepAA_analysis/simulations/scMultiSim/data/phyla3_no_batch_pareto/dataset1/"
# files = "simulation_result_cif_0.5.rds"


# Given a directory with the source simulated data runs MOFA + linear AA on the latent space

run_MOFA <- function(dir){
  
  source("/home/salvatore.milite/data/multideepAA_analysis/simulations/scMultiSim/scripts/utils.R")
  
  library(MOFA2)
  library(dplyr)
  library(Seurat)
  library(Signac)
  library(archetypes)
  library(matrixStats)
  
  
  rds_files <- list.files(dir, pattern = "simulation_result.*[0-9].rds")
  
  for(files in rds_files){
    
    obj <- readRDS(file.path(dir, files))
    metadata <- obj$cell_meta
    n_batches <- metadata$batch %>% unique() %>% length()
    tips <- obj$tree$tip.label
    if(any(grepl("arch",metadata$pop))) {
      n_archetypes <- length(tips)
    } else {
      n_archetypes <- length(tips) + 1
    }
    
    
    obj$counts <- obj$counts[rowVars(obj$counts) >  0 & !is.na(rowVars(obj$counts)) , ]
    obj$atac_counts <- obj$atac_counts[rowVars(obj$atac_counts) > 0 & !is.na(rowVars(obj$atac_counts)), ]
    
    views <-    list(rna = obj$counts   %>% t %>% scale %>% t, 
                     atac = obj$atac_counts %>% t %>% scale %>% t)


    colnames(views$rna) <- paste0("Cell", 1:ncol(views$rna)) 
    colnames(views$atac) <- paste0("Cell", 1:ncol(views$atac)) 
    
    rownames(views$rna) <- paste0("Gene", 1:nrow(views$rna)) 
    rownames(views$atac) <- paste0("Peak", 1:nrow(views$atac)) 
    
    MOFAobject <- create_mofa(views)
 
    data_opts <- get_default_data_options(MOFAobject)
    model_opts <- get_default_model_options(MOFAobject)
    train_opts <- get_default_training_options(MOFAobject)
    
    model_opts$num_factors <- n_archetypes
    
    
    MOFAobject <- prepare_mofa(
      object = MOFAobject,
      data_options = data_opts,
      model_options = model_opts,
      training_options = train_opts
    )
    outfile = file.path(dir, files) %>% gsub(".rds", "_MOFA_real_K.hd5",.)
    mofa_real_K <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)
    outfile_R <- file.path(dir, files) %>% gsub(".rds", "_MOFA_real_K.rds",.)
    mofa_real_K %>% saveRDS(object = ., outfile_R)
    
    train_opts <- get_default_training_options(MOFAobject)
    
    model_opts$num_factors <- n_archetypes - 1
    
    
    MOFAobject <- prepare_mofa(
      object = MOFAobject,
      data_options = data_opts,
      model_options = model_opts,
      training_options = train_opts
    )
    
    outfile = file.path(dir, files) %>% gsub(".rds", "_MOFA_max_K.hd5",.)
    mofa_max_K <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)
    outfile_R <- file.path(dir, files) %>% gsub(".rds", "_MOFA_max_K.rds",.)
    mofa_max_K %>% saveRDS(object = ., outfile_R)
    
    aa_MOFA <- archetypes(data = mofa_max_K@expectations$Z[[1]],
                          k = n_archetypes, verbose = FALSE)
    
    outfile_R <- file.path(dir, files) %>% gsub(".rds", "_MOFA_AA.rds",.)
    aa_MOFA %>% saveRDS(object = ., outfile_R)
    
    
  }
}

# Given a directory with the source simulated data runs intNMF + linear AA on the latent space


run_intNMF <- function(dir){
  
  
  source("/home/salvatore.milite/data/multideepAA_analysis/simulations/scMultiSim/scripts/utils.R")
  
  
  library(dplyr)
  library(IntNMF)
  library(archetypes)
  library(matrixStats)
  
  rds_files <- rds_files <- list.files(dir, pattern = "simulation_result.*[0-9].rds")
  for(files in rds_files){
    
    obj <- readRDS(file.path(dir, files))
    metadata <- obj$cell_meta
    n_batches <- metadata$batch %>% unique() %>% length()
    tips <- obj$tree$tip.label
    if(any(grepl("arch",metadata$pop))) {
      n_archetypes <- length(tips)
    } else {
      n_archetypes <- length(tips) + 1
    }
    
    obj$counts <- obj$counts[rowVars(obj$counts) >  0 & !is.na(rowVars(obj$counts)) , ]
    obj$atac_counts <- obj$atac_counts[rowVars(obj$atac_counts) > 0 & !is.na(rowVars(obj$atac_counts)), ]
    
    views <- list(rna = obj$counts %>%   t  , 
                  atac = obj$atac_counts %>%  t) 
    
    rownames(views$rna) <- paste0("Cell", 1:nrow(views$rna)) 
    rownames(views$atac) <- paste0("Cell", 1:nrow(views$atac)) 
    
    colnames(views$rna) <- paste0("Gene", 1:ncol(views$rna)) 
    colnames(views$atac) <- paste0("Peak", 1:ncol(views$atac)) 
    
    res_NMF_clust <- IntNMF::nmf.mnnals(views, k = n_archetypes, maxiter = 50)    
    res_NMF_max <- IntNMF::nmf.mnnals(views, k = n_archetypes - 1, maxiter = 50)    
    
    
    outfile_R <- file.path(dir, files) %>% gsub(".rds", "_intNMF_real_K.rds",.)
    res_NMF_clust %>% saveRDS(object = ., outfile_R)
    
    outfile_R <- file.path(dir, files) %>% gsub(".rds", "_intNMF_max_K.rds",.)
    res_NMF_max %>% saveRDS(object = ., outfile_R)
    
    aa_intNMF <- archetypes(data = res_NMF_max$W,
                          k = n_archetypes, verbose = FALSE)
    
    outfile_R <- file.path(dir, files) %>% gsub(".rds", "_intNMF_AA.rds",.)
    aa_intNMF %>% saveRDS(object = ., outfile_R)
    
  }
}


run_jive <- function(dir){
  
  
  source("/home/salvatore.milite/data/multideepAA_analysis/simulations/scMultiSim/scripts/utils.R")
  
  
  library(dplyr)
  library(r.jive)
  library(archetypes)
  library(matrixStats)
  
  rds_files <- rds_files <- list.files(dir, pattern = "simulation_result.*[0-9].rds")
  for(files in rds_files){
    obj <- readRDS(file.path(dir, files))
    metadata <- obj$cell_meta
    n_batches <- metadata$batch %>% unique() %>% length()
    
    obj$counts <- obj$counts[rowVars(obj$counts) >  0 & !is.na(rowVars(obj$counts)) , ]
    obj$atac_counts <- obj$atac_counts[rowVars(obj$atac_counts) > 0 & !is.na(rowVars(obj$atac_counts)), ]
    
    
    views <- list(rna = obj$counts %>%  t %>% scale %>% t, 
                  atac = obj$atac_counts   %>% t %>% scale %>% t )
    
    tips <- obj$tree$tip.label
    if(any(grepl("arch",metadata$pop))) {
      n_archetypes <- length(tips)
    } else {
      n_archetypes <- length(tips) + 1
    }
    
    res_jive_clust <- jive(views, rankJ =  n_archetypes, method = "given")
    res_jive_estimated <- jive(views, rankJ =  n_archetypes - 1 , method = "given")
    
    
    outfile_R <- file.path(dir, files) %>% gsub(".rds", "_JIVE_real_K.rds",.)
    res_jive_clust %>% saveRDS(object = ., outfile_R)
    outfile_R <- file.path(dir, files) %>% gsub(".rds", "_JIVE_max_K.rds",.)
    res_jive_estimated %>% saveRDS(object = ., outfile_R)
    pca <- prcomp(res_jive_estimated$joint[[1]] %>% t, scale = TRUE, rank. = res_jive_estimated$rankJ)
    aa_jive <- archetypes(data = pca$x,
                          k = n_archetypes, verbose = FALSE)
    
    outfile_R <- file.path(dir, files) %>% gsub(".rds", "_JIVE_AA.rds",.)
    aa_jive %>% saveRDS(object = ., outfile_R)
    
  }
}

# Given a directory with the source simulated data runs midaa + linear AA on the latent space


run_multideepAA <- function(dir){
  
  
  source("/home/salvatore.milite/data/multideepAA_analysis/simulations/scMultiSim/scripts/utils.R")
  
  
  reticulate::use_condaenv("scdeepaa_sim", required = T)
  library(dplyr)
  library(matrixStats)
  
  midaa <- reticulate::import("midaa")
  
  rds_files <- list.files(dir, pattern = "simulation_result.*[0-9].rds")
  for(files in rds_files){
    obj <- readRDS(file.path(dir, files))
    metadata <- obj$cell_meta
    n_batches <- metadata$batch %>% unique() %>% length()
    
    
    obj$counts <- obj$counts[rowVars(obj$counts) >  0 & !is.na(rowVars(obj$counts)) , ]
    obj$atac_counts <- obj$atac_counts[rowVars(obj$atac_counts) > 0 & !is.na(rowVars(obj$atac_counts)), ]
    
    
    
    views <- list(rna = obj$counts  %>% t %>% scale, 
                  atac = obj$atac_counts  %>% t %>% scale)
    
    
    
    rownames(views$rna) <- paste0("Cell", 1:nrow(views$rna)) 
    rownames(views$atac) <- paste0("Cell", 1:nrow(views$atac)) 
    
    colnames(views$rna) <- paste0("Gene", 1:ncol(views$rna)) 
    colnames(views$atac) <- paste0("Peak", 1:ncol(views$atac)) 
    
    tips <- obj$tree$tip.label
    if(any(grepl("arch",metadata$pop))) {
      n_archetypes <- length(tips)
    } else {
      n_archetypes <- length(tips) + 1
    }
    
    MULT = 64L
    
    res_deepAA_fixZ <- midaa$fit_deepAA(
      list(views$rna, views$atac),
      list(rep(1., dim(views$rna)[1]), rep(1., dim(views$atac)[1])),
      list("G", "G"),
      hidden_dims_dec_common = list(8L * MULT,16L * MULT),
      hidden_dims_dec_last = list(32L * MULT ,64L * MULT),
      hidden_dims_enc_ind = list(64L * MULT, 32L * MULT),
      hidden_dims_enc_common = list(16L * MULT),
      hidden_dims_enc_pre_Z = list(8L * MULT),
      lr = 0.0001,
      gamma_lr = 0.1,
      steps = 2000L,
      narchetypes = as.integer(n_archetypes),
      batch_size = 9000L,
      Z_fix_release_step = as.integer(2000L/3 %>% round),
      Z_fix_norm= 1e-1,
      fix_Z = TRUE,
    )
    inferred <- res_deepAA_fixZ$inferred_quantities$A

    outfile_R <- file.path(dir, files) %>% gsub(".rds", "_deep_Zfix_AA.rds",.)
    res_deepAA_fixZ$inferred_quantities %>% saveRDS(object = ., outfile_R)
    
    res_deepAA_learnZ <- midaa$fit_deepAA(
      list(views$rna, views$atac),
      list(rep(1., dim(views$rna)[1]), rep(1., dim(views$atac)[1])),
      list("G", "G"),
      hidden_dims_dec_common = list(8L * MULT,16L * MULT),
      hidden_dims_dec_last = list(32L * MULT ,64L * MULT),
      hidden_dims_enc_ind = list(64L * MULT, 32L * MULT),
      hidden_dims_enc_common = list(16L * MULT),
      hidden_dims_enc_pre_Z = list(8L * MULT),
      lr = 0.0001,
      gamma_lr = 0.1,
      steps = 1000L,
      narchetypes = as.integer(n_archetypes),
      batch_size = 9000L,
      VAE_steps = 300,
      fix_Z = FALSE
    )
    
    
    outfile_R <- file.path(dir, files) %>% gsub(".rds", "_deep_Zfree_AA.rds",.)
    res_deepAA_learnZ$inferred_quantities %>% saveRDS(object = ., outfile_R)
    
    
    
  }
  return("Done!")
  
}

# Given a directory with the source simulated data runs VAE + linear AA on the latent space


run_multideepAA_auto_only <- function(dir){
  
  
  source("/home/salvatore.milite/data/multideepAA_analysis/simulations/scMultiSim/scripts/utils.R")
  
  
  reticulate::use_condaenv("scdeepaa_sim", required = T)
  library(dplyr)
  library(archetypes)
  library(matrixStats)
  
  midaa <- reticulate::import("midaa")
  
  rds_files <- list.files(dir, pattern = "simulation_result.*[0-9].rds")
  for(files in rds_files){
    obj <- readRDS(file.path(dir, files))
    metadata <- obj$cell_meta
    n_batches <- metadata$batch %>% unique() %>% length()
    
    obj$counts <- obj$counts[rowVars(obj$counts) >  0 & !is.na(rowVars(obj$counts)) , ]
    obj$atac_counts <- obj$atac_counts[rowVars(obj$atac_counts) > 0 & !is.na(rowVars(obj$atac_counts)), ]
    
    
    
    views <- list(rna = obj$counts  %>% t %>% scale, 
                  atac = obj$atac_counts %>% t %>% scale)
    
    rownames(views$rna) <- paste0("Cell", 1:nrow(views$rna)) 
    rownames(views$atac) <- paste0("Cell", 1:nrow(views$atac)) 
    
    colnames(views$rna) <- paste0("Gene", 1:ncol(views$rna)) 
    colnames(views$atac) <- paste0("Peak", 1:ncol(views$atac)) 
    
    tips <- obj$tree$tip.label
    
    if(any(grepl("arch",metadata$pop))) {
      n_archetypes <- length(tips)
    } else {
      n_archetypes <- length(tips) + 1
    }
    
    MULT = 64L
    
    res_auto_narchs <- midaa$fit_deepAA(
      list(views$rna, views$atac),
      list(rep(1., dim(views$rna)[1]), rep(1., dim(views$atac)[1])),
      list("G", "G"),
      hidden_dims_dec_common = list(8L * MULT,16L * MULT),
      hidden_dims_dec_last = list(32L * MULT ,64L * MULT),
      hidden_dims_enc_ind = list(64L * MULT, 32L * MULT),
      hidden_dims_enc_common = list(16L * MULT),
      hidden_dims_enc_pre_Z = list(8L * MULT),
      lr = 0.0001,
      gamma_lr = 0.1,
      steps = 2000L,
      narchetypes = as.integer(n_archetypes) + 1L,
      batch_size = 9000L,
      Z_fix_norm= 1e-2,
      just_VAE = TRUE # We are actually running a VAE
    )
    
    outfile_R <- file.path(dir, files) %>% gsub(".rds", "_multi_realK_AA.rds",.)
    res_auto_narchs$inferred_quantities %>% saveRDS(object = ., outfile_R)
    
    res_auto_hidden_dims <- midaa$fit_deepAA(
      list(views$rna, views$atac),
      list(rep(1., dim(views$rna)[1]), rep(1., dim(views$atac)[1])),
      list("G", "G"),
      hidden_dims_dec_common = list(8L * MULT,16L * MULT),
      hidden_dims_dec_last = list(32L * MULT ,64L * MULT),
      hidden_dims_enc_ind = list(64L * MULT, 32L * MULT),
      hidden_dims_enc_common = list(16L * MULT),
      hidden_dims_enc_pre_Z = list(8L * MULT),
      lr = 0.0001,
      gamma_lr = 0.1,
      steps = 2000L,
      narchetypes = as.integer(n_archetypes),
      batch_size = 9000L,
      fix_Z = FALSE,
      just_VAE = TRUE # We are actually running a VAE
    )
    
    outfile_R <- file.path(dir, files) %>% gsub(".rds", "_multi_maxK_AA.rds",.)
    res_auto_narchs$inferred_quantities %>% saveRDS(object = ., outfile_R)
    
    
    aa_VAE <- archetypes(data = res_auto_hidden_dims$inferred_quantities$Z,
                          k = n_archetypes, verbose = FALSE)
    
    outfile_R <- file.path(dir, files) %>% gsub(".rds", "_multi_AA.rds",.)
    aa_VAE %>% saveRDS(object = ., outfile_R)
    
  }
}

library(easypar)

# Here we use easypar to generate a set of scripts to run as jobs on a SLURM cluster, the package supports
# also PBS and LSF with a similar synthax (https://caravagnalab.github.io/easypar/)

all_dirs <- list.dirs("/home/salvatore.milite/data/multideepAA_analysis/simulations/scMultiSim/data/", recursive = T, full.names = T) %>% 
  grep("dataset[0-9]", ., value = T)  #%>%  grep("phyla3_no_batch_trajectory", ., value = T) #%>% grep("phyla3_no_batch_pareto", ., value = T)

config_cpu <- easypar::default_SBATCH_config()

config_cpu$`--cpus-per-task` <- 4
config_cpu$`--job-name` <- "MOFA_run"
config_cpu$`--mem-per-cpu` <- "16000M"

easypar::run_SLURM(run_MOFA, PARAMS = data.frame(dirs = all_dirs), modules = c("R/4.3.1"),SBATCH_config = config_cpu, output_folder = "../slurm_executables/MOFA_run")

config_cpu$`--cpus-per-task` <- 2

config_cpu$`--job-name` <- "IntNMF_run"


easypar::run_SLURM(run_intNMF, PARAMS = data.frame(dirs = all_dirs),modules = c("R/4.3.1"), SBATCH_config = config_cpu, output_folder = "../slurm_executables/IntNMF_run")

config_cpu$`--cpus-per-task` <- 2

config_cpu$`--job-name` <- "jive_run"

easypar::run_SLURM(run_jive, PARAMS = data.frame(dirs = all_dirs),modules = c("R/4.3.1"),SBATCH_config = config_cpu, output_folder = "../slurm_executables/jive_run")

config_gpu <- easypar::default_SBATCH_config()

config_gpu$`--cpus-per-task` <- 4
config_gpu$`--job-name` <- "deepAA_run"
config_gpu$`--partition` <- "gpuq"
config_gpu$`--time` <- "0:20:00"
config_gpu$`--mem-per-cpu` <- "12000M"
config_gpu$`--gpus` <- 1


easypar::run_SLURM(run_multideepAA, PARAMS = data.frame(dirs = all_dirs),modules = c("R/4.3.1"), SBATCH_config = config_gpu, output_folder = "../slurm_executables/deepAA_run")

config_gpu$`--job-name` <- "VAE_run"

easypar::run_SLURM(run_multideepAA_auto_only, PARAMS = data.frame(dirs = all_dirs),modules = c("R/4.3.1"),  SBATCH_config = config_gpu, output_folder = "../slurm_executables/VAE_run")
