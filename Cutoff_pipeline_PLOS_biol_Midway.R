if (length(commandArgs(trailingOnly=TRUE))==0) {
  args <- c()
} else {
  args <- commandArgs(trailingOnly=TRUE)
}
task <- as.character(args[1]) 
# Can be: edge_weight_distributions | 
# sensitivity_cutoff_selection | 
# sensitivity_cutoff_all_scenarios |
# within_module_diversity


# Functions ---------------------------------------------------------------

source('functions.R')
prep.packages(c('tidyverse','magrittr','sqldf','igraph','data.table','utils'))
scenario_cols <- c('red','orange','blue')
ps_cols <- c('#0A97B7','#B70A97','#97B70A')

gg_labels <- as_labeller(c(`04` = 'Low',
                           `05` = 'Medium',
                           `06` = 'High',
                           `S` = 'Selection',
                           `G` = 'Generalized immunity',
                           `N` = 'Complete neutrality',
                           `Repertoire`='Repertoire',
                           `Module`='Module'))

cutoff_prob_seq = seq(0.25,0.95,0.05)

# Edge weights and cutoffs -----------------------------------------------

if (task=='edge_weight_distributions'){
  print('Getting edge weight distributions...')
  # Selection
  edges_S_04 <- get_edge_disributions(PS = '04',scenario = 'S',exp = '001',1, 0.3)
  edges_S_05 <- get_edge_disributions(PS = '05',scenario = 'S',exp = '001',1, 0.6)
  edges_S_06 <- get_edge_disributions(PS = '06',scenario = 'S',exp = '001',1, 0.85)
  x <- rbind(edges_S_04,edges_S_05,edges_S_06)
  x <- as.tibble(x)
  cutoff_df <- x %>% distinct(PS, cutoff_value, cutoff_prob)
  cutoff_df$label <- paste('PS: ',cutoff_df$PS,'--Cutoff quantile: ',cutoff_df$cutoff_prob,sep='')
  my_labels <- as_labeller(c(`04` = cutoff_df$label[1],
                             `05` = cutoff_df$label[2],
                             `06` = cutoff_df$label[3]))
  png('Results/Edge_weights_distributions_S_R1.png', width = 1920, height = 1080)
  x %>% ggplot(aes(value,fill=PS))+
    geom_density()+
    geom_vline(data=cutoff_df, aes(xintercept = cutoff_value), color='red', size=2)+
    scale_fill_manual(values=ps_cols)+
    facet_wrap(~PS,scales = 'free',labeller = my_labels)+
    theme_bw(base_size=26)
  dev.off()
  # GI
  edges_G_04 <- get_edge_disributions(PS = '04',scenario = 'G',exp = '001',1, 0.3)
  edges_G_05 <- get_edge_disributions(PS = '05',scenario = 'G',exp = '001',1, 0.6)
  edges_G_06 <- get_edge_disributions(PS = '06',scenario = 'G',exp = '001',1, 0.85)
  x <- rbind(edges_G_04,edges_G_05,edges_G_06)
  x <- as.tibble(x)
  cutoff_df <- x %>% distinct(PS, cutoff_value, cutoff_prob)
  cutoff_df$label <- paste('PS: ',cutoff_df$PS,'--Cutoff quantile: ',cutoff_df$cutoff_prob,sep='')
  my_labels <- as_labeller(c(`04` = cutoff_df$label[1],
                             `05` = cutoff_df$label[2],
                             `06` = cutoff_df$label[3]))
  png('Results/Edge_weights_distributions_G_R1.png', width = 1920, height = 1080)
  x %>% ggplot(aes(value,fill=PS))+
    geom_density()+
    geom_vline(data=cutoff_df, aes(xintercept = cutoff_value), color='orange', size=2)+
    scale_fill_manual(values=ps_cols)+
    facet_wrap(~PS,scales = 'free',labeller = my_labels)+
    theme_bw(base_size=26)
  dev.off()
  # Neutral
  edges_N_04 <- get_edge_disributions(PS = '04',scenario = 'N',exp = '001',1, 0.3, get_inter = F)
  edges_N_05 <- get_edge_disributions(PS = '05',scenario = 'N',exp = '001',1, 0.6, get_inter = F)
  edges_N_06 <- get_edge_disributions(PS = '06',scenario = 'N',exp = '001',1, 0.85, get_inter = F)
  x <- rbind(edges_N_04,edges_N_05,edges_N_06)
  x <- as.tibble(x)
  cutoff_df <- x %>% distinct(PS, cutoff_value, cutoff_prob)
  cutoff_df$label <- paste('PS: ',cutoff_df$PS,'--Cutoff quantile: ',cutoff_df$cutoff_prob,sep='')
  my_labels <- as_labeller(c(`04` = cutoff_df$label[1],
                             `05` = cutoff_df$label[2],
                             `06` = cutoff_df$label[3]))
  png('Results/Edge_weights_distributions_N_R1.png', width = 1920, height = 1080)
  x %>% ggplot(aes(value,fill=PS))+
    geom_density()+
    geom_vline(data=cutoff_df, aes(xintercept = cutoff_value), color='blue', size=2)+
    scale_fill_manual(values=ps_cols)+
    facet_wrap(~PS,scales = 'free',labeller = my_labels)+
    theme_bw(base_size=26)
  dev.off()
  
}

# Cutoff sensitivity analysis (SELECTION) --------------------------------------------------------
if (task=='sensitivity_cutoff_selection'){
  print('Getting edge cutoffs for selection only (Boxplots)')
  # Get results
  results_cutoff <- get_results_for_cutoff(cutoff_prob_seq = cutoff_prob_seq, scenario = 'S', run_range = 1:10)
  write_csv(results_cutoff, 'Results/results_cutoff_S.csv')
  
  # Examples for modules
  png('Results/module_examples_S.png', width = 1920, height = 1080)
  results_cutoff %>% 
    filter(run==2) %>%
    distinct(module,layer, PS, cutoff_prob) %>% 
    # filter(scenario=='S') %>% 
    ggplot(aes(x=layer, y=module, group=PS, color=PS))+
    geom_point(size=1)+
    scale_color_manual(values = ps_cols)+
    scale_x_continuous(breaks = seq(0,300,50))+
    labs(y= 'module ID', x='Time (months)', title='Structure example')+
    facet_grid(cutoff_prob~PS, scales='free')+
    mytheme
  dev.off()
  
  # Module and repertoire persistence
  module_persistence <- results_cutoff %>% 
    select(PS, run, cutoff_prob, layer, module) %>% 
    group_by(PS,run,cutoff_prob,module) %>% 
    summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
    mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
    mutate(type='Module') %>% 
    rename(id=module) %>% 
    mutate(id=as.character(id))
  
  strain_persistence <- results_cutoff %>% 
    select(PS, run, cutoff_prob, layer, strain_cluster) %>% 
    group_by(PS,run,cutoff_prob,strain_cluster) %>% 
    summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
    mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
    mutate(type='Repertoire') %>% 
    rename(id=strain_cluster) %>% 
    mutate(id=as.character(id))
  
  
  persistence_df <- module_persistence %>% bind_rows(strain_persistence) 
  write_csv(persistence_df, 'Results/persistence_df_cutoff_S.csv')
  
  png('Results/persistence_boxplots.png', width = 1920, height = 1080)
  persistence_df %>% 
    ggplot(aes(x=cutoff_prob, y=persistence, group=cutoff_prob, fill=PS))+
    geom_boxplot(outlier.shape = NA)+
    stat_summary(fun.y=mean, color="red", geom="point", 
                 shape=18, size=3,show.legend = FALSE)+
    facet_grid(PS~type, scales='free', labeller = gg_labels)+
    scale_x_continuous(breaks = cutoff_prob_seq)+
    labs(x='Cut off', y=' Persistence')+
    scale_fill_manual(values=ps_cols)+mytheme
  dev.off()
  
  png('Results/relative_persistence_boxplots.png', width = 1920, height = 1080)
  persistence_df %>% 
    ggplot(aes(x=cutoff_prob, y=relative_persistence, group=cutoff_prob, fill=PS))+
    geom_boxplot(outlier.shape = NA)+
    stat_summary(fun.y=mean, color="red", geom="point", 
                 shape=18, size=3,show.legend = FALSE)+  
    facet_grid(PS~type, scales='free', labeller = gg_labels)+
    scale_x_continuous(breaks = cutoff_prob_seq)+
    labs(x='Cut off', y=' Relative persistence')+
    scale_fill_manual(values=ps_cols)+mytheme
  dev.off()
  
  # Modules per layer
  png('Results/modules_per_layer_boxplots.png', width = 1920, height = 1080)
  results_cutoff %>% 
    group_by(PS, scenario, run, cutoff_prob, layer) %>% 
    summarise(modules_per_layer=length(unique(module))) %>% 
    ggplot(aes(x=cutoff_prob, y=modules_per_layer, group=cutoff_prob, fill=PS))+
    geom_boxplot(outlier.shape = NA)+
    stat_summary(fun.y=mean, color="red", geom="point", 
                 shape=18, size=3,show.legend = FALSE)+
    facet_grid(~PS, scales='free', labeller = gg_labels)+
    scale_x_continuous(breaks = cutoff_prob_seq)+
    labs(x='Cut off', y=' Modules per layer')+
    scale_fill_manual(values=ps_cols)+mytheme
  dev.off()
  
  # Reps per module
  png('Results/Repertoires_per_module_boxplots.png', width = 1920, height = 1080)
  results_cutoff %>% 
    group_by(PS, scenario, run, cutoff_prob, module) %>% 
    summarise(repertoires_per_module=length(unique(strain_cluster))) %>% 
    # group_by(PS, scenario, run, cutoff_prob) %>% summarise(repertoires_per_module=mean(repertoires_per_module)) %>%
    ggplot(aes(x=cutoff_prob, y=repertoires_per_module, group=cutoff_prob, fill=PS))+
    geom_boxplot(outlier.size=0)+
    stat_summary(fun.y=mean, color="red", geom="point", 
                 shape=18, size=3,show.legend = FALSE)+
    facet_grid(~PS, scales='free', labeller = gg_labels)+
    scale_x_continuous(breaks = cutoff_prob_seq)+
    labs(x='Cut off', y='Repertoires per module')+
    scale_fill_manual(values=ps_cols)+mytheme
  dev.off()
  
  
  # quantiles_95 <- function(x) {
  #   r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  #   names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  #   r
  # }
  # results_cutoff %>% 
  #   group_by(PS, scenario, run, cutoff_prob, module) %>% 
  #   summarise(repertoires_per_module=length(unique(strain_cluster))) %>% 
  #   # group_by(PS, scenario, run, cutoff_prob) %>% summarise(repertoires_per_module=mean(repertoires_per_module)) %>%
  #   ggplot(aes(x=cutoff_prob, y=repertoires_per_module, group=cutoff_prob, fill=PS))+
  #   stat_summary(fun.data = quantiles_95, geom="boxplot")+
  #   stat_summary(fun.y=mean, color="red", geom="point", 
  #                shape=18, size=3,show.legend = FALSE)+
  #   facet_grid(~PS, scales='free', labeller = ps_labels)+
  #   scale_x_continuous(breaks = cutoff_prob_seq)+
  #   labs(x='Cut off', y='Repertoires per module')+
  #   scale_fill_manual(values=ps_cols)+mytheme
  
}

# Cutoff sensitivity analysis (all scenarios) --------------------------------------------
if (task=='sensitivity_cutoff_all_scenarios'){
  print('Getting edge distributions for all scenarios...')
  
  results_cutoff_N <- get_results_for_cutoff(cutoff_prob_seq = cutoff_prob_seq, scenario = 'N', run_range = 1:10)
  write_csv(results_cutoff_N, 'Results/results_cutoff_N.csv')
  module_persistence <- results_cutoff_N %>% 
    select(PS, run, cutoff_prob, layer, module) %>% 
    group_by(PS,run,cutoff_prob,module) %>% 
    summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
    mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
    mutate(type='Module') %>% 
    rename(id=module) %>% 
    mutate(id=as.character(id))
  strain_persistence <- results_cutoff_N %>% 
    select(PS, run, cutoff_prob, layer, strain_cluster) %>% 
    group_by(PS,run,cutoff_prob,strain_cluster) %>% 
    summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
    mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
    mutate(type='Repertoire') %>% 
    rename(id=strain_cluster) %>% 
    mutate(id=as.character(id))
  persistence_df_N <- module_persistence %>% bind_rows(strain_persistence) 
  write_csv(persistence_df_N, 'Results/persistence_df_cutoff_N.csv')
  
  
  results_cutoff_G <- get_results_for_cutoff(cutoff_prob_seq = cutoff_prob_seq, scenario = 'G', run_range = 1:10)
  write_csv(results_cutoff_G, 'Results/results_cutoff_G.csv')
  module_persistence <- results_cutoff_G %>% 
    select(PS, run, cutoff_prob, layer, module) %>% 
    group_by(PS,run,cutoff_prob,module) %>% 
    summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
    mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
    mutate(type='Module') %>% 
    rename(id=module) %>% 
    mutate(id=as.character(id))
  strain_persistence <- results_cutoff_G %>% 
    select(PS, run, cutoff_prob, layer, strain_cluster) %>% 
    group_by(PS,run,cutoff_prob,strain_cluster) %>% 
    summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
    mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
    mutate(type='Repertoire') %>% 
    rename(id=strain_cluster) %>% 
    mutate(id=as.character(id))
  persistence_df_G <- module_persistence %>% bind_rows(strain_persistence) 
  write_csv(persistence_df_G, 'Results/persistence_df_cutoff_G.csv')
  
  
  # Join the scenarios to one data frame
  results_cutoff_S <- read_csv('Results/results_cutoff_S.csv', col_types = 'iicccicccid')
  results_cutoff <- rbind(results_cutoff_S,results_cutoff_N,results_cutoff_G)
  results_cutoff %<>% mutate(scenario=factor(scenario, levels=c('S','G','N')))
  
  persistence_df_S <- read_csv('Results/persistence_df_cutoff_S.csv', col_types = 'cidciiidc')
  persistence_df <- rbind(persistence_df_S,persistence_df_G,persistence_df_N)
  persistence_df %<>% mutate(scenario=factor(scenario, levels=c('S','G','N')))
  
  
  
  png('Results/module_persistence_scenarios.png', width = 1920, height = 1080)
  persistence_df %>% 
    filter(type=='Module') %>% 
    ggplot(aes(x=cutoff_prob, color=scenario))+
    geom_point(aes(y=mean_persistence))+geom_line(aes(y=mean_persistence))+
    geom_point(aes(y=median_persistence), shape=15)+geom_line(aes(y=median_persistence), linetype='dashed')+
    facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
    scale_x_continuous(breaks = cutoff_prob_seqs)+
    labs(x='Cut off', y=' Mean or median MODULE persistence')+
    scale_color_manual(values=scenario_cols)+mytheme
  dev.off()
  
  png('Results/repertoire_persistence_scenarios.png', width = 1920, height = 1080)
  persistence_df %>% 
    filter(type=='Repertoire') %>% 
    ggplot(aes(x=cutoff_prob, color=scenario))+
    geom_point(aes(y=mean_persistence))+geom_line(aes(y=mean_persistence))+
    geom_point(aes(y=median_persistence), shape=15)+geom_line(aes(y=median_persistence), linetype='dashed')+
    facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
    scale_x_continuous(breaks = cutoff_prob_seq)+
    labs(x='Cut off', y=' Mean or median REPERTOIRE persistence')+
    scale_color_manual(values=scenario_cols)+mytheme
  dev.off()
  
  
  png('Results/module_relative_persistence_scenarios.png', width = 1920, height = 1080)
  persistence_df %>% 
    filter(type=='Module') %>% 
    ggplot(aes(x=cutoff_prob, color=scenario))+
    geom_point(aes(y=mean_relative_persistence))+geom_line(aes(y=mean_relative_persistence))+
    geom_point(aes(y=median_relative_persistence), shape=15)+geom_line(aes(y=median_relative_persistence), linetype='dashed')+
    facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
    scale_x_continuous(breaks = cutoff_prob_seq)+
    labs(x='Cut off', y=' Mean or median MODULE relative persistence')+
    scale_color_manual(values=scenario_cols)+mytheme
  dev.off()
  
  png('Results/repertoire_relative_persistence_scenarios.png', width = 1920, height = 1080)
  persistence_df %>% 
    filter(type=='Repertoire') %>% 
    ggplot(aes(x=cutoff_prob, color=scenario))+
    geom_point(aes(y=mean_relative_persistence))+geom_line(aes(y=mean_relative_persistence))+
    geom_point(aes(y=median_relative_persistence), shape=15)+geom_line(aes(y=median_relative_persistence), linetype='dashed')+
    facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
    scale_x_continuous(breaks = cutoff_prob_seq)+
    labs(x='Cut off', y=' Mean or median REPERTOIRE relative persistence')+
    scale_color_manual(values=scenario_cols)+mytheme
  dev.off()
  
  # Modules per layer
  png('Results/modules_per_layer_scenarios.png', width = 1920, height = 1080)
  results_cutoff %>% 
    group_by(PS, scenario, run, cutoff_prob, layer) %>% 
    summarise(modules_per_layer=length(unique(module))) %>% 
    group_by(PS, scenario, cutoff_prob) %>%
    summarise(mean_modules_per_layer=mean(modules_per_layer),
              median_modules_per_layer=median(modules_per_layer)) %>%
    ggplot(aes(x=cutoff_prob, color=scenario))+
    geom_point(aes(y=mean_modules_per_layer))+geom_line(aes(y=mean_modules_per_layer))+
    geom_point(aes(y=median_modules_per_layer), shape=15)+geom_line(aes(y=median_modules_per_layer), linetype='dashed')+
    facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
    scale_x_continuous(breaks = cutoff_prob_seq)+
    labs(x='Cut off', y=' Mean or median modules per layer')+
    scale_color_manual(values=scenario_cols)+mytheme
  dev.off()
  
  # Reps per module
  png('Results/repertoires_per_module_scenarios.png', width = 1920, height = 1080)
  results_cutoff %>%
    group_by(PS, scenario, run, cutoff_prob, module) %>%
    summarise(repertoires_per_module=length(unique(strain_cluster))) %>%
    group_by(PS, scenario, cutoff_prob) %>%
    summarise(mean_repertoires_per_module=mean(repertoires_per_module),
              median_repertoires_per_layer=median(repertoires_per_module)) %>%
    ggplot(aes(x=cutoff_prob, color=scenario))+
    geom_point(aes(y=mean_repertoires_per_module))+geom_line(aes(y=mean_repertoires_per_module))+
    geom_point(aes(y=median_repertoires_per_layer), shape=15)+geom_line(aes(y=median_repertoires_per_layer), linetype='dashed')+
    facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
    scale_x_continuous(breaks = cutoff_prob_seq)+
    labs(x='Cut off', y=' Mean or median Repertoires per module')+
    scale_color_manual(values=scenario_cols)+mytheme
  dev.off()
  
}

# Within-module diversity -------------------------------------------------
if (task=='within_module_diversity'){
  print ('Calculating within-module temporal diversity...')
  
  for (scenario in c('S','N','G')){
    design_cutoff <- expand.grid(PS=sprintf('%0.2d', 4:6),
                                 scenario=scenario, 
                                 exp='001',
                                 run_range=1:10,
                                 cutoff_prob=cutoff_prob_seq,
                                 stringsAsFactors = F)
    statistic_results <- c()
    for (i in 1:nrow(design_cutoff)){
      PS <- design_cutoff[i,1]
      scenario <- design_cutoff[i,2]
      exp <- design_cutoff[i,3]
      run <- design_cutoff[i,4]
      cutoff_prob <- design_cutoff[i,5]
      x <- calculate_module_diversity(PS, scenario, exp, run, cutoff_prob)
      statistic_results <- rbind(statistic_results, x)
      # print(paste(object.size(statistic_results), nrow(statistic_results)))
    }
    write_csv(statistic_results, paste('Results/temporal_diversity_cutoff_',scenario,'.csv',sep=''))
    
    # statistic_results <- read_csv(paste('Results/statistic_results_',scenario,'.csv',sep=''))
    
    print(paste('Plotting statistic for scenario',scenario))
    png(paste('Results/temporal_diversity_cutoff_',scenario,'.png',sep=''), width = 1920, height = 1080)
    statistic_results %>% 
      ggplot(aes(x=cutoff_prob, y=statistic, group=cutoff_prob, fill=PS))+
      geom_boxplot(outlier.size=0)+
      stat_summary(fun.y=mean, color=scenario_cols[which(c('S','N','G')==scenario)], geom="point", 
                   shape=18, size=3,show.legend = FALSE)+
      facet_grid(~PS, scales='free', labeller = gg_labels)+
      scale_x_continuous(breaks = cutoff_prob_seq)+
      labs(x='Cut off', y='Temporal diversity')+
      scale_fill_manual(values=ps_cols)+mytheme
    dev.off()
  }
  
  # Now plot for all scenarios
  s <- read_csv('Results/temporal_diversity_cutoff_S.csv', col_types = 'ccidiiiiddd')
  n <- read_csv('Results/temporal_diversity_cutoff_N.csv', col_types = 'ccidiiiiddd')
  g <- read_csv('Results/temporal_diversity_cutoff_G.csv', col_types = 'ccidiiiiddd')
  
  statistic_results <- rbind(s,n,g)
  
  x <- statistic_results %>% 
    mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% 
    group_by(scenario, PS, cutoff_prob) %>% 
    summarise(mean_D=mean(D),
              median_D=median(D),
              sd_D=sd(D),
              mean_statistic=mean(statistic),
              median_statistic=median(statistic),
              sd_statistic=sd(statistic),
              n=length(statistic)) %>% 
    mutate(ymin_D_sd=mean_D-sd_D, 
           ymax_D_sd=mean_D+sd_D,
           ymin_D_ci=mean_D-qnorm(0.975)*sd_D/sqrt(n),
           ymax_D_ci=mean_D+qnorm(0.975)*sd_D/sqrt(n),
           ymin_statistic_sd=mean_statistic-sd_statistic, 
           ymax_statistic_sd=mean_statistic+sd_statistic,
           ymin_statistic_ci=mean_statistic-qnorm(0.975)*sd_statistic/sqrt(n),
           ymax_statistic_ci=mean_statistic+qnorm(0.975)*sd_statistic/sqrt(n))
  
  png('Results/D_scenarios.png', width = 1920, height = 1080)
  x %>% 
    ggplot(aes(x=cutoff_prob, color=scenario))+
    geom_point(aes(y=mean_D))+geom_line(aes(y=mean_D))+
    geom_point(aes(y=median_D), shape=15)+geom_line(aes(y=median_D), linetype='dashed')+
    geom_errorbar(aes(ymin=ymin_D_ci,ymax=ymax_D_ci),alpha=0.6, width=0.01, size=1)+
    facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
    scale_x_continuous(breaks = cutoff_prob_seq)+
    labs(x='Cut off', y=' Mean or median D')+
    scale_color_manual(values=scenario_cols)+mytheme
  dev.off()
  
  png('Results/temporal_diversity_cutoff_scenarios.png', width = 1920, height = 1080)
  x %>% 
    ggplot(aes(x=cutoff_prob, color=scenario))+
    geom_point(aes(y=mean_statistic))+geom_line(aes(y=mean_statistic))+
    geom_point(aes(y=median_statistic), shape=15)+geom_line(aes(y=median_statistic), linetype='dashed')+
    geom_errorbar(aes(ymin=ymin_statistic_ci,ymax=ymax_statistic_ci),alpha=0.6, width=0.01, size=1)+
    facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
    scale_x_continuous(breaks = cutoff_prob_seq)+
    labs(x='Cut off', y=' Mean or median Temporal Diversity')+
    scale_color_manual(values=scenario_cols)+mytheme
  dev.off()
  
}
