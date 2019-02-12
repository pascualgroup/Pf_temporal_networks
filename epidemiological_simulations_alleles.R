# Initialize --------------------------------------------------------------
if (length(commandArgs(trailingOnly=TRUE))==0) {
  args <- c('06','S','001',0.85,300,12,10,2)
} else {
  print('Taking arguments from command line.')
  args <- commandArgs(trailingOnly=TRUE)
}
PS <- as.character(args[1])
scenario <- as.character(args[2])
exp <- as.character(args[3])
run <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
cutoff_prob <- as.numeric(args[4])
numLayers <- as.numeric(args[5]) # This is to limit the number of layers. When running the real model this should be at the maximum value. number of layers will be (MaxTime-18000)/window width, whre 18000 is the burnin time of the model and window width is in days (typically 30)
time_interval <- as.numeric(args[6])
n_hosts <- as.numeric(args[7]) # Number of naive hosts to infect
n_samples <- as.numeric(args[8]) # Number of random starting point layers within each module

# print('Arguments passed:')
# cat('experiment: ');cat(experiment);cat('\n')
# cat('run: ');cat(run);cat('\n')
# cat('numLayers: ');cat(numLayers);cat('\n')
# cat('cutoffPercentile: ');cat(cutoffPercentile);cat('\n')
# cat('time interval: ');cat(time_interval);cat('\n')
# cat('# hosts: ');cat(n_hosts);cat('\n')
# cat('# samples: ');cat(n_samples);cat('\n')

# if (Sys.info()[4]=='ee-pascual-dell01'){source('/home/shai/Documents/malaria_temporal_networks/mtn_functions.R')}

source('functions.R')
prep.packages(c("tidyverse","data.table",'magrittr','sqldf'))

base_name <- paste('PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,sep='')
print(base_name)
print('Loading modules and strain compositon...')
modules <- infomap_readTreeFile(PS, scenario, exp, run, cutoff_prob, paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,sep=''))

sampled_strains <- modules$sampled_strains
sampled_alleles <- modules$sampled_alleles
sampled_alleles$allele_locus <- paste(sampled_alleles$allele, sampled_alleles$locus,sep='_')
sampled_alleles %<>% select(gene_id, allele_locus) %>% arrange(gene_id,allele_locus)
sampled_strains %<>% left_join(sampled_alleles)
sampled_strains$strain_id <- as.character(sampled_strains$strain_id )
all(modules$modules$strain_id %in% sampled_strains$strain_id)

print('Merging the strain composition with the module data frame...')
rep_module_layer <- left_join(modules$modules, sampled_strains, by='strain_id') %>% distinct(strain_id,module,layer)
  
duration_new_allele <- 360/120 # 360 is the naive duration of infection and 120 is the number of alleles in a repertoire, both taken from the experiment parameters run in the agent-based model

# How many infections per layer should be performed? That would depend on the 
# EIR parameter from the stochastic ABM which is on the units of
# infectious bites/person/month. So infections_per_layer=length_of_layer*EIR
tmp <- get_data(parameter_space = PS, scenario = scenario, experiment = exp, run = run, cutoff_prob = cutoff_prob, use_sqlite = T, tables_to_get = 'summary_general')[[1]]
infections_per_layer <- ceiling(mean(tmp$EIR)) # need to use ceiling because for low diversity/transmission EIR is < 1.

print('Finished initializing')

# Functions ---------------------------------------------------------------

sample_min_distance <- function(L, d, N){
  # Function originally programmed in MatLab taken from here: https://stackoverflow.com/questions/31971344/generating-random-sequence-with-minimum-distance-between-elements-matlab
  #L is length of interval
  #d is minimum distance
  #N is number of points
  
  E=L-(N-1)*d;  #excess space for points
  
  P=L-d
  while (max(P)>=L-d){
    #generate N+1 random values; 
    Ro=runif(N+1)   # random vector
    #normalize so that the extra space is consumed
    #extra value is the amount of extra space "unused"
    Rn=E*Ro[1:N]/sum(Ro); #normalize
    
    #spacing of points
    S=d*matrix(rep(1,N),nrow = N,ncol=1)+Rn;  
    
    #location of points, adjusted to "start" at 0
    P=round(cumsum(S)-1)
  }
  if(any(diff(P)<d)){
    print('Failed generating sequence. returning NULL')
    return(NULL)
  } else {
    return(P)
  }
}


# A function to follow event queues and produce infections
simulate_infections <- function(event){
  repertoires_seen <- c()
  alleles_seen <- c()
  infection_id <- 0
  infection_queue <- data.frame(module=NULL, layer=NULL, infection_id=NULL, alleles_seen=NULL, repertoires_seen=NULL, duration=NULL)
  for (i in 1:nrow(event)){
    infection_id <- infection_id+1
    rep <- event$strain_id[i]
    rep_alleles <- unique(subset(sampled_strains,strain_id==rep)$allele_locus) # Alleles in the repertoire
    new_alleles <- sum(!rep_alleles%in%alleles_seen) # Alleles in that repertorie which are new to the host
    duration <- new_alleles*duration_new_allele
    # update history
    alleles_seen <- unique(c(alleles_seen, rep_alleles))
    repertoires_seen <- unique(c(repertoires_seen, rep))
    # Add to infection history queue
    infection_queue[infection_id,'layer'] <- event[infection_id,'layer']
    infection_queue[infection_id,'module'] <- event[infection_id,'module']
    infection_queue[infection_id,'infection_id'] <- infection_id
    infection_queue[infection_id,'alleles_seen'] <- length(alleles_seen)
    infection_queue[infection_id,'repertoires_seen'] <- length(repertoires_seen)
    infection_queue[infection_id,'duration'] <- duration
  }
  return(infection_queue)
}

event_status <- function(events){
  print(paste('Total events:',nrow(events)))
  print(paste('Mean infections per layer:',mean(table(events$layer))))
  print(paste('Number of distinct modules:',length(unique(events$module))))
  print(paste('Number of distinct repertoires:',length(unique(events$strain_id))))
}


build_event_queue_within_modules <- function(){
  # This selects n_hosts modules, each with time_interval layers
  x <- rep_module_layer %>% group_by(module,layer) %>% summarise(nreps=length(strain_id))
  # Select the first 12 layers in each module. The trick when using top_n is to give the layer as the nverse weight to give priority to lower numbers
  event_queue <- x %>% group_by(module) %>% top_n(n=time_interval, wt=1/layer) %>% arrange(module,layer) 
  
  # Remain only with those modules which have at least time_interval layers
  y <- event_queue %>% group_by(module) %>% summarise(n=n()) %>% filter(n>=time_interval) # calculate the number of layers each module is present
  event_queue <- inner_join(event_queue,y)
  # Limit to the number of modules
  if(length(unique(event_queue$module))<n_hosts){
    message('Number of modules smaller than n_hosts; leaving all of them')
    print('Number of modules smaller than n_hosts; leaving all of them')
  }
  if(length(unique(event_queue$module))>n_hosts){
    m <- sort(sample(unique(event_queue$module), n_hosts, F))
    event_queue <- event_queue %>% filter(module%in%m)
  }
  
  # Which repertoires appear in these module-layer combinations?
  x <- rep_module_layer %>% filter(module%in%event_queue$module & layer%in%event_queue$layer) %>% group_by(module,layer)%>% arrange(module,layer) 
  event_queue <- inner_join(x, event_queue) %>% arrange(module,layer) 
  # Select infections_per_layer. Notice that there is replacement because this
  # is within module and some modules do not have enough repertories within a
  # layer. That is the whole idea in the within-module case...
  event_queue <- event_queue %>% group_by(module, layer) %>% sample_n(infections_per_layer, T)
  
  verify <- event_queue %>% group_by(module) %>% summarise(count=length(unique(layer)))
  if(!all(verify$count==time_interval)){warning('Something went wrong; apparently not all module-layer combinations have the correct time_interval')}
  
  return(event_queue)
}


build_event_queue_between_modules <- function(){
  
  # In this case we need to select blocks of time_interval infections with
  # at least infections_per_layer in each layer
  # Then select 1 repertoire from each module in that layer
  # And limit the number of repertoires to that determined by biting rate
  
  starting_points <- sort(sample_min_distance(L=numLayers, d=time_interval+1, N=n_hosts)) # sample_min_distance makes sure there will be no overlap in layer sequences
  layers <- unlist(lapply(starting_points, function(x) seq(from=x,to=x+time_interval-1,by=1)))
  # select 1 repertoire from infections_per_layer modules in each layer. Select
  # with replacement in case there are less than infections_per_layer modules in
  # a layer.
  x <- subset(rep_module_layer, layer%in%layers)
  event_queue <- NULL
  for (l in layers){
    y <- subset(x, layer==l)
    if (length(unique(y$module))<infections_per_layer){
      z <- y %>% group_by(layer) %>% sample_n(infections_per_layer, replace = T)
    } else {
      z <- y %>% group_by(module,layer) %>% sample_n(1) # Select 1 repertoire from each module
      z <- z %>% group_by(layer) %>% sample_n(infections_per_layer, replace = F) # Select infections_per_layer repertoires
    }  
    event_queue <- rbind(event_queue, z)
  }
  return(event_queue)
}


build_event_queue_random <- function(){
  
  # In this case is similar to that of between-modules but we randomly
  # select infections_per_layer from each layer, regardless of modlules
  
  starting_points <- sort(sample_min_distance(L=numLayers, d=time_interval+1, N=n_hosts)) # sample_min_distance makes sure there will be no overlap in layer sequences
  layers <- unlist(lapply(starting_points, function(x) seq(from=x,to=x+time_interval-1,by=1)))
  # select infections_per_layer repertoires in each layer, with replacement in
  # case there are less than infections_per_layer modules in a layer.
  x <- subset(rep_module_layer, layer%in%layers)
  
  event_queue <- x %>% group_by(layer) %>% sample_n(infections_per_layer, replace = T) %>% arrange(layer,module)
  
  return(event_queue)
}


# Run simulatinos ---------------------------------------------------------

# results_n_samples <- NULL
for (s in 1:n_samples){
  # Within-module simulations 
  # In the case of within-module infections each host is actually a module because
  # we follow a host for 12 layers within a module
  print(paste(base_name,' | sample ',s,' | within',sep=''))
  events_within <- build_event_queue_within_modules() # build event queue
  event_status(events_within)
  infection_history_within <- NULL
  for (e in unique(events_within$module)){
    print(paste(Sys.time(),' | ',base_name,' | sample ',s,' | event ',which(unique(events_within$module)==e),' | within',sep=''))
    x <- subset(events_within, module==e)
    y <- simulate_infections(x)
    infection_history_within <- bind_rows(infection_history_within, y)
  }
  # infection_history_within %>% group_by(infection_id) %>% summarise(d=mean(duration)) %>% ggplot(aes(infection_id, d))+geom_point()+geom_line()
  
  
  # Between-module simulations
  # Each host is a sequence of 10 consecutive layers
  events_between <- build_event_queue_between_modules() # build event queue
  event_status(events_between)
  ## Split the queue to single events (hosts)
  d <- unique(events_between$layer)
  x <- seq_along(d)
  host_events <-  split(d, ceiling(x/time_interval))
  infection_history_between <- NULL
  for (e in 1:n_hosts){
    print(paste(Sys.time(),' | ',base_name,' | sample ',s,' | host ',e,' | between',sep=''))
    # print(e)
    x <- subset(events_between, layer%in%host_events[[e]])
    y <- simulate_infections(x)
    infection_history_between <- bind_rows(infection_history_between, y)
  }
  # infection_history_between %>% group_by(infection_id) %>% summarise(d=mean(duration)) %>% ggplot(aes(infection_id, d))+geom_point()+geom_line()
  
  
  # Random-module simulations
  # Each host is a sequence of 10 consecutive layers
  events_random <- build_event_queue_random() # build event queue
  event_status(events_random)
  ## Split the queue to single events (hosts)
  d <- unique(events_random$layer)
  x <- seq_along(d)
  host_events <-  split(d, ceiling(x/time_interval))
  infection_history_random <- NULL
  for (e in 1:n_hosts){
    print(paste(Sys.time(),' | ',base_name,' | sample ',s,' | host ',e,' | random',sep=''))
    x <- subset(events_random, layer%in%host_events[[e]])
    y <- simulate_infections(x)
    infection_history_random <- bind_rows(infection_history_random, y)
  }
  # infection_history_random %>% group_by(infection_id) %>% summarise(d=mean(duration)) %>% ggplot(aes(infection_id, d))+geom_point()+geom_line()
  
  # Finalize
  ## Join the results
  infection_history_between$case <- 'B'
  infection_history_within$case <- 'W'
  infection_history_random$case <- 'R'
  results <- bind_rows(infection_history_between,infection_history_within,infection_history_random)
  results$sample <- s
  results$PS <- PS
  results$scenario <- scenario
  results$exp <- exp
  results$run <- run
  results$cutoff_prob <- cutoff_prob

    ## Write results to file
  write.table(results, paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/epi_alleles/',base_name,'_epi_T',time_interval,'_H',n_hosts,'_sample_',s,'.csv',sep=''), sep=',')
  
  # Write the event queues
  events_within$case <- 'W'
  events_between$case <- 'B'
  events_random$case <- 'R'
  events <- bind_rows(events_within[,c("strain_id","module","layer","case")],events_between,events_random)
  events$sample <- s
  events$PS <- PS
  events$scenario <- scenario
  events$exp <- exp
  events$run <- run
  events$cutoff_prob <- cutoff_prob
  write.table(events, paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/epi_alleles/',base_name,'_epi_T',time_interval,'_H',n_hosts,'_sample_',s,'_EVENTS.csv',sep=''), sep=',')
  
  # results_n_samples <- rbind(results_n_samples, results)
}
# 
# results_n_samples_S <- results_n_samples
# 
# as_tibble(results) %>%
#   bind_rows(results_n_samples_S) %>%
#   group_by(scenario, case, infection_id) %>%
#   summarise(mean_d=mean(duration)) %>%
#   ggplot(aes(infection_id, mean_d, color=scenario))+
#   geom_point()+
#   geom_line()+
#   facet_wrap(~case)
