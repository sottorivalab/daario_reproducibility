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
    
    rep = as.integer(files %>% gsub( ".*_rep_", "", .) %>% gsub( ".rds", "", ., fixed = TRUE) )
    n_cells = as.integer(files %>% gsub( ".*_ncells_", "", .) %>% gsub( "_rep_.*.rds", "", .) )
    
    obj$counts <- obj$counts[rowVars(obj$counts) >  0 & !is.na(rowVars(obj$counts)) , ]
    obj$atac_counts <- obj$atac_counts[rowVars(obj$atac_counts) > 0 & !is.na(rowVars(obj$atac_counts)), ]
    
    
    
    views <- list(rna = obj$counts  %>% t %>% scale, 
                  atac = obj$atac_counts %>% t %>% scale)
    
    rept <- n_cells / nrow(views$rna)
    
    views$rna <- replicate(expr = views$rna,rept, simplify = F) %>% do.call(rbind,.)
    views$atac <- replicate(expr = views$atac,rept, simplify = F) %>% do.call(rbind,.)
    
    rownames(views$rna) <- paste0("Cell", 1:nrow(views$rna)) 
    rownames(views$atac) <- paste0("Cell", 1:nrow(views$atac)) 
    
    colnames(views$rna) <- paste0("Gene", 1:ncol(views$rna)) 
    colnames(views$atac) <- paste0("Peak", 1:ncol(views$atac)) 
    
    tips <- obj$tree$tip.label
    n_archetypes <- length(tips)
    
    MULT = 8L
    
    prova <- midaa$fit_DAARIO(
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
      steps = 2L,
      narchetypes = n_archetypes,
      batch_size = 4096L,
      
    )
    
    start_time <- Sys.time()
    
    res_deepAA_learnZ <- midaa$fit_DAARIO(
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
      steps = 500L,
      narchetypes = n_archetypes,
      batch_size = 4096L,

    )
    end_time <- Sys.time()
    
    name_to_substitute <- "_timing_results.rds"
    outfile_R <- file.path(dir, files) %>% gsub(".rds", name_to_substitute,.)
    data.frame(
      N_cells = n_cells,
      rep = rep, 
      time = end_time - start_time 
    ) %>% saveRDS(object = ., outfile_R)
    
    
    
  }
  return("Done!")
  
}

library(dplyr)

config_gpu <- easypar::default_SBATCH_config()

config_gpu$`--cpus-per-task` <- 4
config_gpu$`--job-name` <- "deepAA_run_timing"
config_gpu$`--partition` <- "gpuq"
config_gpu$`--time` <- "4:00:00"
config_gpu$`--mem-per-cpu` <- "12000M"
config_gpu$`--gpus` <- 1



all_dirs <- list.dirs("/home/salvatore.milite/data/multideepAA_analysis/simulations/scMultiSim/data_timing/", recursive = T, full.names = T) %>% 
  grep("dataset[0-9]", ., value = T) 



easypar::run_SLURM(run_multideepAA, PARAMS = data.frame( dirs = all_dirs ),
                   modules = c("R/4.3.1"), per_task = 1,
                   SBATCH_config = config_gpu, output_folder = "../slurm_executables_timing/deepAA_run")