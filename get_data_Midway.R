# Initialize --------------------------------------------------------------
source('functions.R')
# prep.packages(c('sqldf','tidyverse','magrittr','igraph','data.table'))
library(sqldf, quietly = T, warn.conflicts = F)
library(tidyverse, quietly = T, warn.conflicts = F)
library(magrittr, quietly = T, warn.conflicts = F)
library(igraph, quietly = T, warn.conflicts = F)
library(data.table, quietly = T, warn.conflicts = F)

if (length(commandArgs(trailingOnly=TRUE))==0) {
  args <- c('550','S','001',0.85, '118,126,138,142,154,162',3)
  # args <- c('550','S','001',0.85, '10:100')
} else {
  args <- commandArgs(trailingOnly=TRUE)
}
job_ps <- as.character(args[1]) 
job_scenario <- as.character(args[2]) 
job_exp <- as.character(args[3]) 
# job_run <- 1
job_run <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
cutoff_prob <- as.numeric(args[4])
layers <- as.character(args[5]) # Layers can be in a 1:300 format for a sequence or in a '118,126,138,142,154,162' format for separate layers.
task <- as.character(args[6]) # Can be: make_networks | prepare_infomap | read_infomap_results | temporal_diversity | module_Fst
modularity_exp <-  as.numeric(args[7]) # 0: not comparing to empirical; 1: genes+MOI=1; 2: alleles+MOI=1; 3: alleles, all repertoires; 4: genes, all repertoires
write_edge_weights <- as.character(args[8]) # Should write the edge weights before cutoff?
  
if (str_detect(layers,'\\:')){
  layers <- seq(str_split(layers,'\\:')[[1]][1],str_split(layers,'\\:')[[1]][2],1)
} else {
  layers <- as.integer(str_split(layers,",")[[1]])  
}

base_name <- paste('PS',job_ps,'_',job_scenario,'_E',job_exp,'_R',job_run,'_',cutoff_prob,sep='')

print(base_name)
print(layers)
print(task)
print(modularity_exp)
if (modularity_exp!=0){print('Comparing to empirical data')}

# Make networks -----------------------------------------------------------
make_network <- function(unit_for_edges, repertoires_to_sample, write_to_files, write_edge_weights){
  
  # Data from sqlite
  data <- get_data(parameter_space = job_ps, scenario = job_scenario, experiment = job_exp, run = job_run, host_age_structure = T, use_sqlite = T)
  write_csv(data$summary_general, paste(base_name,'_summary_general.csv',sep=''))
  write_csv(data$sampled_infections, paste(base_name,'_sampled_infections.csv',sep=''))
  
  # Network objects
  if (modularity_exp!=0){
    cutoff_value <- NULL
  }
  if(modularity_exp==0){
    # If experiment is not control then take the cutoff value from the control. This
    # requires the control to be run first.
    if (job_exp=='001' | !file.exists(paste('PS',job_ps,'_',job_scenario,'_E001_R',job_run,'_',cutoff_prob,'_network_info.csv',sep=''))){
      cutoff_value <- NULL
    } else {
      x <- readLines(paste('PS',job_ps,'_',job_scenario,'_E001_R',job_run,'_',cutoff_prob,'_network_info.csv',sep=''))
      cutoff_value <- x[6]
    }
  }
  
  
  print('Building network...')
  message('Building network...')
  network <- createTemporalNetwork(ps = job_ps,
                                   scenario = job_scenario,
                                   exp = job_exp, 
                                   run = job_run,
                                   cutoff_prob = cutoff_prob, 
                                   cutoff_value = cutoff_value,
                                   layers_to_include = layers,
                                   sampled_infections = data$sampled_infections,
                                   unit_for_edges = unit_for_edges,
                                   repertoires_to_sample = repertoires_to_sample,
                                   write_to_files=write_to_files)
  
  # edges <- c(network$intralayer_edges_no_cutoff, network$interlayer_edges_no_cutoff)
  # as.tibble(edges) %>% ggplot(aes(value))+geom_density()+geom_vline(xintercept = network$cutoff_value)
  
  sink(paste(base_name,'_network_info.csv',sep=''), append = F)
  cat(network$ps);cat('\n')
  cat(network$scenario);cat('\n')
  cat(network$experiment);cat('\n')
  cat(network$run);cat('\n')
  cat(network$cutoff_prob);cat('\n')
  cat(network$cutoff_value);cat('\n')
  sink.reset()
  write_csv(network$layer_summary, paste(base_name,'_layer_summary.csv',sep=''))
  
  if(write_edge_weights){
    write_csv(as.tibble(network$intralayer_edges_no_cutoff), paste(base_name,'_intralayer_no_cutoff.csv',sep=''))
    write_csv(as.tibble(network$interlayer_edges_no_cutoff), paste(base_name,'_interlayer_no_cutoff.csv',sep=''))
  }
  
  print('Finished building network...')
  message('Finished building network...')
  return(network)
}

if (task=='make_networks'){
  make_network(unit_for_edges = 'alleles', repertoires_to_sample = NULL, write_to_files = F, write_edge_weights = write_edge_weights)
}



# Strain persistence with usearch -----------------------------------------

if (task=='repertoire_persistence'){
  
  print('Getting information on genes and alleles...')
  # Get the gene list for the repertoires
  if (on_Midway()){
    db <- dbConnect(SQLite(), dbname = paste('/scratch/midway2/pilosofs/PLOS_Biol/sqlite/','PS',job_ps,'_',job_scenario,'_E',job_exp,'_R',job_run,'.sqlite',sep=''))
  } else {
    db <- dbConnect(SQLite(), dbname = paste('/media/Data/PLOS_Biol/sqlite_',job_scenario,'/','PS',job_ps,'_',job_scenario,'_E',job_exp,'_R',job_run,'.sqlite',sep=''))
  }
  
  sampled_strains <- as.tibble(dbGetQuery(db, 'SELECT id,gene_id FROM sampled_strains'))
  names(sampled_strains)[1] <- 'strain_id'
  sampled_alleles <- as.tibble(dbGetQuery(db, 'SELECT * FROM sampled_alleles'))
  sampled_infections <- as.tibble(dbGetQuery(db, 'SELECT * FROM sampled_infections'))
  dbDisconnect(db)
  
  # Find unique repertoires using usearch. This needs to be done because the ABM
  # does not keep track on strain composition to give strains unique names.
  print('Finding unique repertoires using usearch...')
  
  # Produce allele combinations in repertoires
  sampled_alleles$allele_locus <- paste(sampled_alleles$allele, sampled_alleles$locus,sep='_')
  sampled_alleles %<>% select(gene_id, allele_locus) %>% arrange(gene_id,allele_locus)
  sampled_strains %<>% left_join(sampled_alleles) %>% distinct(strain_id,allele_locus)
  
  # Make a file for usearch. Need to use letters because it does not accept 0 and 1
  strain_allele_mat <- xtabs(~strain_id+allele_locus, sampled_strains)
  strain_allele_mat[strain_allele_mat==0] <- 'A'
  strain_allele_mat[strain_allele_mat==1] <- 'G'
  f <- paste(base_name,'_strain_sequences_without_modules.fasta',sep='')
  # f <- paste('/media/Data/PLOS_Biol/parameter_files/',base_name,'_strain_sequences_without_modules.fasta',sep='')
  sink(f, append = F)
  for (i in 1:nrow(strain_allele_mat)){
    cat('>');cat(rownames(strain_allele_mat)[i]);cat('\n')
    cat(paste(strain_allele_mat[i,],collapse=''));cat('\n')
  }
  sink.reset()
  # rep_clusters_file <- paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/',job_ps,'_',job_scenario,'/',base_name,'_unique_repertoires_without_modules.txt',sep='')
  rep_clusters_file <- paste(base_name,'_unique_repertoires_without_modules.txt',sep='')
  system(paste("/scratch/midway2/pilosofs/PLOS_Biol/usearch10.0.240_i86linux32 -fastx_uniques ",f," -relabel rep -tabbedout ",rep_clusters_file,sep='')) # Find unique sequences
  # system(paste("/media/Data/PLOS_Biol/parameter_files/usearch10.0.240_i86linux32 -fastx_uniques ",f," -relabel rep -tabbedout ",rep_clusters_file,sep='')) # Find unique sequences
  
  print('Getting results from usearch...')
  unique_repertoires <- read_delim(rep_clusters_file,
                                   delim='\t', 
                                   col_names = c('strain_id','strain_cluster','cluster','member_id','repertoires_in_cluster','first_strain_id'))
  
  unique_repertoires$strain_id <- as.character(unique_repertoires$strain_id)
  sampled_strains$strain_id <- as.character(sampled_strains$strain_id)
  sampled_strains %<>% left_join(unique_repertoires) %>% distinct(strain_id,strain_cluster)
  
  print('Calculate repertorie persistence...')
  sampled_infections %<>% distinct(strain_id, time)
  sampled_infections$strain_id <- as.character(sampled_infections$strain_id)
  sampled_strains %<>% left_join(sampled_infections)
  
  layer_df <- tibble(time=sort(unique(sampled_strains$time)))
  layer_df$layer <- 1:nrow(layer_df)
  
  sampled_strains %<>% left_join(layer_df)
  
  sampled_strains %<>% 
    group_by(strain_cluster) %>% 
    summarise(birth_layer=min(layer), death_layer=max(layer)) %>% 
    mutate(persistence=death_layer-birth_layer+1)
  
  write_csv(sampled_strains, paste(base_name,'_repertoire_persistence_without_modules.txt',sep=''))
}

# Build infomap objects ---------------------------------------------------


if (task=='prepare_infomap'){

  if (modularity_exp==0){
    network <- make_network(unit_for_edges = 'alleles', repertoires_to_sample = NULL, write_to_files = F, write_edge_weights = write_edge_weights) # First make the network
  }
  if (modularity_exp==1){
    network <- make_network(unit_for_edges = 'genes', repertoires_to_sample = c(90,68,69,52,115,44), write_to_files = F, write_edge_weights = write_edge_weights) # First make the network
    # network <- make_network(unit_for_edges = 'genes', repertoires_to_sample = c(90,66,65,55), write_to_files = T, write_edge_weights = T) # First make the network
  }
  if (modularity_exp==2){
    # network <- make_network(unit_for_edges = 'alleles', repertoires_to_sample = c(98,68,69,52), write_to_files = T, write_edge_weights = T) # First make the network
    network <- make_network(unit_for_edges = 'alleles', repertoires_to_sample = c(90,68,69,52,115,44), write_to_files = F, write_edge_weights = write_edge_weights) # First make the network
  }
  if (modularity_exp==3){
    network <- make_network(unit_for_edges = 'alleles', repertoires_to_sample = NULL, write_to_files = F, write_edge_weights = write_edge_weights) # First make the network
  }
  if (modularity_exp==4){
    network <- make_network(unit_for_edges = 'genes', repertoires_to_sample = NULL, write_to_files = F, write_edge_weights = write_edge_weights) # First make the network
  } 
  
  print('Building Infomap objects...')
  message('Building Infomap objects...')
  # Should interlayer edges be rescaled? Only if working with the 6 layers of the
  # interventions in the PS of the interventions.
  if (modularity_exp!=0){
    infomap <- build_infomap_objects(network_object = network, 
                                     write_to_infomap_file = T, 
                                     infomap_file_name = paste(base_name,'_Infomap_multilayer.txt',sep=''), 
                                     return_objects = T,
                                     repertoire_survival_prob = NULL,
                                     rescale_by_survival_prob=F
                                    )
  } else {
    infomap <- build_infomap_objects(network_object = network, 
                                     write_to_infomap_file = T, 
                                     infomap_file_name = paste(base_name,'_Infomap_multilayer.txt',sep=''), 
                                     return_objects = T,
                                     rescale_by_survival_prob=F
                                     )
  }
  write_csv(infomap$nodeList, paste(base_name,'_node_list.csv',sep=''))
} # End task 'prepare_infomap'


# Read Infomap results ----------------------------------------------------

if (task=='read_infomap_results'){
  print('Getting results from Infomap...')
  x <- infomap_readTreeFile(PS = job_ps,scenario = job_scenario,exp = job_exp,run=job_run,cutoff_prob = cutoff_prob)
  
  # Find unique repertoires using usearch. This needs to be done because the ABM
  # does not keep track on strain composition to give strains unique names.
  print('Finding unique repertoires using usearch...')
  modules <- x$modules
  sampled_strains <- x$sampled_strains
  sampled_alleles <- x$sampled_alleles
  
  print(paste('Identical strains in module and strain data frames?',setequal(sampled_strains$strain_id, modules$strain_id)))
  print(paste('Are all the strains in the module data frame contained in the strain data frame?',all(modules$strain_id%in%sampled_strains$strain_id)))
  
  # Produce allele combinations in repertoires
  sampled_alleles$allele_locus <- paste(sampled_alleles$allele, sampled_alleles$locus,sep='_')
  sampled_alleles %<>% select(gene_id, allele_locus) %>% arrange(gene_id,allele_locus)
  sampled_strains %<>% left_join(sampled_alleles) %>% distinct(strain_id,allele_locus)
  
  # Only consider the strains that appear in the modules. That is particularly
  # necessary for cases where not all 300 layers are analyzed for modularity
  # (like IRS)
  sampled_strains %<>% filter(strain_id%in%modules$strain_id)
  
  # Make a file for usearch. Need to use letters because it does not accept 0 and 1
  strain_allele_mat <- xtabs(~strain_id+allele_locus, sampled_strains)
  strain_allele_mat[strain_allele_mat==0] <- 'A'
  strain_allele_mat[strain_allele_mat==1] <- 'G'
  f <- paste(base_name,'_strain_sequences.fasta',sep='')
  # f <- paste('/media/Data/PLOS_Biol/parameter_files/',base_name,'_strain_sequences.fasta',sep='')
  sink(f, append = F)
  for (i in 1:nrow(strain_allele_mat)){
    cat('>');cat(rownames(strain_allele_mat)[i]);cat('\n')
    cat(paste(strain_allele_mat[i,],collapse=''));cat('\n')
  }
  sink.reset()
  # rep_clusters_file <- paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/',job_ps,'_',job_scenario,'/',base_name,'_unique_repertoires.txt',sep='')
  rep_clusters_file <- paste(base_name,'_unique_repertoires.txt',sep='')
  system(paste("/scratch/midway2/pilosofs/PLOS_Biol/usearch10.0.240_i86linux32 -fastx_uniques ",f," -relabel rep -tabbedout ",rep_clusters_file,sep='')) # Find unique sequences
  # system(paste("/media/Data/PLOS_Biol/parameter_files/usearch10.0.240_i86linux32 -fastx_uniques ",f," -relabel rep -tabbedout ",rep_clusters_file,sep='')) # Find unique sequences
  
  print('Getting results from usearch...')
  unique_repertoires <- read_delim(rep_clusters_file,
                                   delim='\t', 
                                   col_names = c('strain_id','strain_cluster','cluster','member_id','repertoires_in_cluster','first_strain_id'))
  
  unique_repertoires$strain_id <- as.character(unique_repertoires$strain_id)
  sampled_strains$strain_id <- as.character(sampled_strains$strain_id)
  modules$strain_id <- as.character(modules$strain_id)
  
  sampled_strains %<>% left_join(unique_repertoires) %>% select(strain_id,allele_locus,strain_cluster)
  modules %<>% left_join(unique_repertoires) %>% 
    mutate(strain_cluster=parse_number(strain_cluster)) %>% 
    select(layer, module, strain_cluster, strain_id, strain_unique, nodeID, PS, scenario, exp, run, cutoff_prob, -first_strain_id, -repertoires_in_cluster, -member_id, -cluster)
    
  
  print('Writing Infomap results with unique repertoires to files...')
  write_csv(sampled_strains, paste(base_name,'_sampled_strains.csv',sep=''))
  write_csv(sampled_alleles, paste(base_name,'_sampled_alleles.csv',sep=''))
  
  # # When not analyzing all layers (like in IRS), then before writing the module
  # # data frame, change the layer index numbers to the actual layers.
  # if (layers!=1:300){
  #   layers_idx <- tibble(layer=1:length(layers),real_layer=layers)
  #   modules %<>% left_join(layers_idx) %>% select(-layer) %>% rename(layer=real_layer) %>% select(layer, module:cutoff_prob)
  #   # If want to keep the layer_idx do that. But it will have an additional column, which will create problems with compartibility for functions that use modules.csv
  #   # modules %<>% left_join(layers_idx) %>% rename(layer_idx=layer) %>% rename(layer=real_layer) %>% select(layer, module:cutoff_prob, layer_idx)
  # }
  write_csv(modules, paste(base_name,'_modules.csv',sep=''))
  print('Done!')
} # End task read_infomap_results


# Calcualte diversity -----------------------------------------------------


if (task=='temporal_diversity'){
  print('Calculating module temporal diversity')
  # Module temporal diversity
  x <- calculate_module_diversity(job_ps,job_scenario,job_exp,job_run,cutoff_prob)
  write_csv(x, paste(base_name,'_temporal_diversity.csv',sep=''))
}

if (task=='module_Fst'){
  print('Calculating module Fst')
  #mFst
  x <- calculate_mFst(job_ps,job_scenario,job_exp,job_run,cutoff_prob, maxModule = 75)
  write_csv(x, paste(base_name,'_mFst.csv',sep=''))
}
