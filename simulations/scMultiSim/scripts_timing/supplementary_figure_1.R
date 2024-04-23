# Load all data
library(tidyverse)
library(ggbeeswarm)

extrafont::loadfonts(quiet = T)

all_fls <- list.files("../data_timing/", recursive = T, pattern = "_timing_results.rds", full.names = T)

dfs <- lapply(all_fls, readRDS) %>% do.call(rbind,.)


dfs %>%  ggplot(aes(x = as.numeric(N_cells), y = as.numeric(time), group = N_cells)) +
  geom_boxplot() +
  geom_beeswarm() +
  
  #geom_line(linewidth = 1., alpha = 0.6, color = "grey10") +
  #stat_summary(fun.y = "median", geom = "point", size = 2.5, fun.max = "max", fun.min = "min", shape = 1) +
  theme_classic(base_size = 12, base_family = "Arial") + scale_y_log10() + scale_x_log10() + ylab("seconds") + xlab("# cells") +
  Seurat::RotatedAxis()
