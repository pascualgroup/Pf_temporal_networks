# Initialize --------------------------------------------------------------
source('functions.R')
# prep.packages(c('sqldf','tidyverse','magrittr','igraph','data.table'))
library(sqldf, quietly = T, warn.conflicts = F)
library(tidyverse, quietly = T, warn.conflicts = F)
library(magrittr, quietly = T, warn.conflicts = F)
library(igraph, quietly = T, warn.conflicts = F)
library(data.table, quietly = T, warn.conflicts = F)
library(utils, quietly = T, warn.conflicts = F)


if (length(commandArgs(trailingOnly=TRUE))==0) {
  args <- c('S','001')
} else {
  args <- commandArgs(trailingOnly=TRUE)
}
PS <- Sys.getenv('SLURM_ARRAY_TASK_ID')
PS <- str_pad(PS, 2, 'left', '0')
scenario <- as.character(args[1]) 
exp <- as.character(args[2]) 

print(paste(PS,scenario,exp,sep=' | '))

experiments <- as.tibble(expand.grid(PS=PS,
                           scenario=scenario, 
                           exp=exp,
                           run=1:50,
                           stringsAsFactors = F))
cutoff_df <- tibble(PS=sprintf('%0.2d', c(4:6,18)),
                    cutoff_prob=c(0.3,0.6,0.85,0.85)
                    )
cutoff_df$folder <- paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/',cutoff_df$PS,'_',scenario,'/',sep='')
print(cutoff_df)

experiments <- experiments %>% left_join(cutoff_df)

print(experiments)

module_results <- c()
for (i in 1:nrow(experiments)){
  run <- experiments$run[i]
  exp <- experiments$exp[i]
  cutoff_prob <- experiments$cutoff_prob[i]
  folder <- experiments$folder[i]
  x <- get_modularity_results(PS = PS,scenario = scenario,exp = exp,run = run,cutoff_prob = cutoff_prob,folder = folder)
  module_results <- rbind(module_results,x)
}
module_results$strain_cluster <- as.integer(module_results$strain_cluster)
write_csv(module_results, paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/PS',PS,'_',scenario,'_E',exp,'_module_results.csv',sep=''))

# Relative persistence
module_persistence <- module_results %>% 
  select(scenario, PS, run, cutoff_prob, layer, module) %>% 
  group_by(scenario, PS,run,cutoff_prob,module) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Module') %>% 
  rename(id=module) %>% mutate(id=as.character(id))

strain_persistence <- module_results %>% 
  select(scenario, PS, run, cutoff_prob, layer, strain_cluster) %>% 
  group_by(scenario, PS,run,cutoff_prob,strain_cluster) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Repertoire') %>% 
  rename(id=strain_cluster) %>% mutate(id=as.character(id))

persistence_df <- module_persistence %>% bind_rows(strain_persistence)
write_csv(persistence_df, paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/PS',PS,'_',scenario,'_E',exp,'_persistence_df.csv',sep=''))

# Temporal diversity
temporal_diversity <- c()
for (i in 1:nrow(experiments)){
  run <- experiments$run[i]
  cutoff_prob <- experiments$cutoff_prob[i]
  exp <- experiments$exp[i]
  folder <- experiments$folder[i]
  x <- get_temporal_diversity(PS = PS,scenario = scenario,exp = exp,run = run,cutoff_prob = cutoff_prob,folder = folder)
  temporal_diversity <- rbind(temporal_diversity,x)
}
temporal_diversity <- as.tibble(temporal_diversity)
write_csv(temporal_diversity, paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/PS',PS,'_',scenario,'_E',exp,'_temporal_diversity.csv',sep=''))

# mFst
mFst <- c()
for (i in 1:nrow(experiments)){
  run <- experiments$run[i]
  exp <- experiments$exp[i]
  cutoff_prob <- experiments$cutoff_prob[i]
  folder <- experiments$folder[i]
  x <- get_mFst(PS = PS,scenario = scenario,exp = exp,run = run,cutoff_prob = cutoff_prob,folder = folder)
  mFst <- rbind(mFst,x)
}
mFst <- as.tibble(mFst) 
write_csv(mFst, paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/PS',PS,'_',scenario,'_E',exp,'_mFst.csv',sep=''))
