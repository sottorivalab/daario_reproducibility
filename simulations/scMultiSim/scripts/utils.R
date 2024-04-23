
library(scales)
library(dplyr)

groupwise_means <- function(group_name, df, mat, round = F) {
  group_cells <- df$cell_id[df$pop == group_name]
  group_data <- mat[, group_cells, drop=FALSE]
  group_mean <- rowMeans(group_data, na.rm = TRUE)
  if(round) group_mean <- round(group_mean)
  return(group_mean)
}


# A set of functions that just load the A and B matrices for the different input methods

load_all_MOFA <- function(dir = "/group/sottoriva/salvatore.milite/multideepAA_analysis/simulations/scMultiSim/data/phyla3_no_batch_pareto/") {
  
  files_right_K <- list.files(dir, pattern = "MOFA", recursive = T, full.names = T) %>% grep("_real_K.rds", .,value = T)
  files_AA <- list.files(dir, pattern = "MOFA", recursive = T, full.names = T) %>% grep("_AA.rds", .,value = T)
  
  data_right_K <- lapply(files_right_K, FUN = function(x) {
    obj <- readRDS(x)
    return(obj@expectations$Z[[1]])
  })
  names(data_right_K) <- files_right_K
  
  data_AA <- lapply(files_AA, FUN = function(x) {
    obj <- readRDS(x)
    return(obj[c("alphas", "betas", "archetypes")])
  })
  names(data_AA) <- files_AA
  
  return(list(AA = data_AA, real_K = data_right_K))
}


load_all_deepAA <- function(dir = "/group/sottoriva/salvatore.milite/multideepAA_analysis/simulations/scMultiSim/data/phyla3_no_batch_pareto/") {
  
  files_AA_freeZ <- list.files(dir, pattern = "_deep", recursive = T, full.names = T) %>% grep("Zfree", .,value = T)
  files_AA_fixZ <- list.files(dir, pattern = "_deep", recursive = T, full.names = T) %>% grep("Zfix", .,value = T)
  
  data_AA_freeZ <- lapply(files_AA_freeZ, FUN = function(x) {
    obj <- readRDS(x)
    return(obj[c("A", "B", "Z")])
  })
  names(data_AA_freeZ) <- files_AA_freeZ
  
  data_AA_fixZ <- lapply(files_AA_fixZ, FUN = function(x) {
    obj <- readRDS(x)
    return(obj[c("A", "B", "Z")])
  })
  names(data_AA_fixZ) <- files_AA_fixZ
  
  return(list(AA_Zfix = data_AA_fixZ, AA_Zfree = data_AA_freeZ))
}


load_all_VAE <- function(dir = "/group/sottoriva/salvatore.milite/multideepAA_analysis/simulations/scMultiSim/data/phyla3_no_batch_pareto/") {
  
  files_right_K <- list.files(dir, pattern = "_multi", recursive = T, full.names = T) %>% grep("_realK_AA", .,value = T)
  files_AA <- list.files(dir, pattern = "_multi", recursive = T, full.names = T) %>% grep("multi_AA.rds", .,value = T)
  
  data_right_K <- lapply(files_right_K, FUN = function(x) {
    obj <- readRDS(x)
    return(obj$Z)
  })
  names(data_right_K) <- files_right_K
  
  data_AA <- lapply(files_AA, FUN = function(x) {
    obj <- readRDS(x)
    return(obj[c("alphas", "betas", "archetypes")])
  })
  names(data_AA) <- files_AA
  
  return(list(AA = data_AA, real_K = data_right_K))
}




load_all_NMF <- function(dir = "/group/sottoriva/salvatore.milite/multideepAA_analysis/simulations/scMultiSim/data/phyla3_no_batch_pareto/") {
  
  files_right_K <- list.files(dir, pattern = "intNMF", recursive = T, full.names = T) %>% grep("_real_K", .,value = T)
  files_AA <- list.files(dir, pattern = "intNMF", recursive = T, full.names = T) %>% grep("_AA.rds", .,value = T)
  
  data_right_K <- lapply(files_right_K, FUN = function(x) {
    obj <- readRDS(x)
    return(obj$W)
  })
  names(data_right_K) <- files_right_K
  
  data_AA <- lapply(files_AA, FUN = function(x) {
    obj <- readRDS(x)
    return(obj[c("alphas", "betas", "archetypes")])
  })
  names(data_AA) <- files_AA
  
  return(list(AA = data_AA, real_K = data_right_K))
}

load_all_jive <- function(dir = "/group/sottoriva/salvatore.milite/multideepAA_analysis/simulations/scMultiSim/data/phyla3_no_batch_pareto/") {
  
  files_right_K <- list.files(dir, pattern = "_JIVE", recursive = T, full.names = T) %>% grep("_real_K", .,value = T)
  files_AA <- list.files(dir, pattern = "_JIVE", recursive = T, full.names = T) %>% grep("_AA.rds", .,value = T)
  
  data_right_K <- lapply(files_right_K, FUN = function(x) {
    obj <- readRDS(x)
    return(prcomp(obj$joint[[1]], scale = TRUE, rank. = obj$rankJ)$rotation)
  })
  names(data_right_K) <- files_right_K
  
  data_AA <- lapply(files_AA, FUN = function(x) {
    obj <- readRDS(x)
    return(obj[c("alphas", "betas", "archetypes")])
  })
  names(data_AA) <- files_AA
  
  return(list(AA = data_AA, real_K = data_right_K))
}


load_all_input <- function(dir = "/group/sottoriva/salvatore.milite/multideepAA_analysis/simulations/scMultiSim/data/phyla3_no_batch_pareto/") {
  
  all_files <- list.files(dir, pattern = "[0-9].rds$", recursive = T, full.names = T) 

  data_all <- lapply(all_files, FUN = function(x) {
    obj <- readRDS(x)
    return(obj)
  })
  names(data_all) <- all_files
  
  return(data = data_all)
}



compute_euclidean_distance <- function(col1, col2) {
  # Compute the Euclidean distance
  sqrt(sum((col1 - col2)^2, na.rm = T))
}

rename_cols_based_on_distance <- function(df1, df2) {
  # Ensure data is numeric
  df1 <- data.matrix(df1)
  df2 <- data.matrix(df2)
  
  # Compute the Euclidean distance between all column pairs
  distance_matrix <- pairwise_euclidean(df1, df2)
  
  # Initialize new names vector
  new_names <- vector("character", ncol(df2))
  
  for(col_index in 1:ncol(df2)) {
    
    col_index <- which.min(colMins(distance_matrix, na.rm = T))
    

    # Extract the distances related to the current df2 column
    current_distances <- distance_matrix[, col_index]
    
    # Get the name of the df1 column with the smallest distance
    min_index <- which.min(current_distances)
    new_names[col_index] <- colnames(df1)[min_index]
    
    # Remove already chosen names from consideration by setting their distance artificially high
    already_chosen <- which(colnames(df1) %in% new_names)
    #current_distances[already_chosen] <- Inf
    distance_matrix[already_chosen,] <- Inf
    distance_matrix[,col_index] <- Inf
  }
  
  # Rename columns
  colnames(df2) <- new_names
  
  return(df2)
}


compute_euclidean_distance_df <- function(df1, df2, aggregation_fun = mean) {
  # Compute the Euclidean distance
  df1 <- df1[, sort(colnames(df1))]
  df2 <- df2[, sort(colnames(df1))]
  to_rt <- sapply(1:nrow(df1), function(i) compute_euclidean_distance(df1[i,],  df2[i,]))
  return(aggregation_fun(to_rt))
}

compute_ARI_df <- function(df1, df2, aggregation_fun = mean) {
  # Compute the Euclidean distance
  df1 <- df1[, sort(colnames(df1))]
  df2 <- df2[, sort(colnames(df1))]
  
  clusters_df1 <- apply(df1, 1, which.max) %>% as.character()
  clusters_df2 <- apply(df2, 1, which.max) %>% as.character()
  
  return(aricode::ARI(clusters_df1, clusters_df2))
}

pairwise_euclidean <- function(matrix1, matrix2) {
  # Ensure the inputs are correct
  if (ncol(matrix1) == 0 || ncol(matrix2) == 0) {
    stop("Both matrices must have at least one column.")
  }
  
  if (nrow(matrix1) != nrow(matrix2)) {
    stop("Both matrices must have the same number of rows.")
  }
  
  # Get the number of columns in each matrix
  num_cols_matrix1 <- ncol(matrix1)
  num_cols_matrix2 <- ncol(matrix2)
  
  # Initialize a matrix to hold the results
  distances <- matrix(0, nrow = num_cols_matrix1, ncol = num_cols_matrix2)
  
  # Calculate Euclidean distances
  for (i in 1:num_cols_matrix1) {
    for (j in 1:num_cols_matrix2) {
      distances[i, j] <- sqrt(sum((matrix1[, i] - matrix2[, j])^2, na.rm = T))
    }
  }
  
  rownames(distances) <- colnames(matrix1)
  colnames(distances) <- colnames(matrix2)
  
  # Return the matrix of distances
  distances
}


calculate_norm_logCPM <- function(count_matrix, log = F) {
  
  # Step 1: Calculate the library size for each cell
  library_sizes <- colSums(count_matrix)
  
  # Step 2: Calculate the scaling factor for each cell (as per million counts)
  scaling_factors <- library_sizes / 1e3
  
  # Step 3: Divide each count by the scaling factor for its cell (normalize)
  # This accounts for differences in sequencing depth between cells
  norm_count_matrix <- sweep(count_matrix, 2, scaling_factors, "/")

  # Step 4: Add a pseudocount to avoid undefined log-values for zero counts
  pseudocount <- 0.5
  norm_count_matrix_pseudo <- norm_count_matrix + pseudocount
  
  # Step 5: Apply log transformation (natural logarithm)
  if(log){
    log_norm_count_matrix <- log(norm_count_matrix_pseudo)
  } else{
    log_norm_count_matrix <- round(norm_count_matrix)
  }


  return(log_norm_count_matrix)
}


get_extreme_states_from_trajectory_old <- function(df, annot, tree) {
  annot_extrema <- annot %>%  separate(pop, into = c("from", "pop"), sep = "_") %>% 
    filter(pop %in% tree$tip.label) %>% select(-from)
  extreme_points <- annot_extrema %>% group_by(pop) %>% top_frac(.05, wt = depth) 
  colnames(df) <- paste0("cell", 1:ncol(df))
  first_points <- annot %>% top_frac(.05, wt = -depth) 
  first_points$pop <- "init"
  df_filt <- df[,c(extreme_points$cell_id, first_points$cell_id)]
  df_filt_split <- lapply(c(extreme_points$pop %>% unique(), "init"),
                          groupwise_means,df =  rbind(extreme_points,
                                                      first_points), mat = df_filt)  %>% do.call(rbind,.) %>% t
  colnames(df_filt_split) <- paste0("arc", 1:ncol(df_filt_split))
  return(df_filt_split)
}

get_extreme_states_from_trajectory <- function(df, annot, tree) {
  annot_extrema <- annot %>%  separate(pop, into = c("from", "pop"), sep = "_") %>% 
    filter(pop %in% tree$tip.label) %>% select(-from)
  extreme_points <- annot_extrema %>% group_by(pop) %>% top_frac(.15, wt = depth) 
  colnames(df) <- paste0("cell", 1:ncol(df))
  df_filt <- df[,c(extreme_points$cell_id)]
  df_filt_split <- lapply(c(extreme_points$pop %>% unique()),
                          groupwise_means,df =  extreme_points,
                                                       mat = df_filt)  %>% do.call(rbind,.) %>% t
  colnames(df_filt_split) <- paste0("arc", 1:ncol(df_filt_split))
  return(df_filt_split)
}


get_extreme_states_from_AA <- function(A, counts, meta) {
  
  colnames(A) <- paste0("arc", 1:ncol(A))
  meta_arch <- cbind(meta, A) %>% select(-pop, -depth)
  extreme_points <- list()
  for(i in colnames(A)) {
    extreme_points[[i]] <-  meta_arch %>% select(cell_id, !!i)  %>% 
                          top_frac(.15) %>% pull(cell_id) 
  }
  
  rownames(counts) <- paste0("cell", 1:nrow(counts))
  df_filt_split <- lapply(extreme_points,
                          function(x) {
                            colMeans(counts[x,])
                          })  %>% do.call(rbind,.) %>% t
  colnames(df_filt_split) <- paste0("arc", 1:ncol(df_filt_split))
  return(df_filt_split)
  
  
}
