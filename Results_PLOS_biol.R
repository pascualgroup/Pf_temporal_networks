## @knitr Initialize

# Functions ---------------------------------------------------------------

source('~/Documents/malaria_interventions/functions.R')
prep.packages(c('tidyverse','magrittr','sqldf','igraph','data.table','googlesheets','utils','cowplot','grid','gridExtra'), verbose = F)

# Initialize important variables ------------------------------------------
setwd('/media/Data/PLOS_Biol/')

monitored_variables <- c('prevalence', 'meanMOI','n_circulating_strains', 'n_circulating_genes', 'n_alleles', 'n_total_bites')

scenario_cols <- c('red','blue','orange') # Order is: S, N, G
ps_cols <- c('#0A97B7','#B70A97','#97B70A','#873600')
gg_labels <- as_labeller(c(`04` = 'Low',
                           `05` = 'Medium',
                           `06` = 'High',
                           `18` = 'Seasonal (High)',
                           `S` = 'Selection',
                           `G` = 'Generalized immunity',
                           `N` = 'Complete neutrality',
                           `prevalence` = 'Prevalence',
                           `meanMOI` = 'MOI (mean)',
                           `n_circulating_strains` = '# repertoires',
                           `n_circulating_genes` = '# genes',
                           `n_alleles` = '# alleles',
                           `n_total_bites` = '# bites',
                           `Module` = 'Module',
                           `Repertoire` = 'Repertoire'))

all_experiments <- expand.grid(PS=sprintf('%0.2d', c(4:6,18)),
                           scenario=c('S','N','G'), 
                           exp='001',
                           run=1:50,
                           # cutoff_prob=seq(0.3,0.95,0.05),
                           stringsAsFactors = F)
cutoff_df <- tibble(PS=sprintf('%0.2d', c(4:6,18)),cutoff_prob=c(0.3,0.6,0.85,0.85))
all_experiments <- left_join(all_experiments,cutoff_df)


# files_list <- list.files('/media/Data/PLOS_Biol/Results/cutoff_to_use')
# files <- tibble(file=files_list,
#                 PS = sapply(str_split(files_list,'_'),function (x) parse_number(x[1])),
#                 scenario=sapply(str_split(files_list,'_'),function (x) x[2]),
#                 experiment=sapply(str_split(files_list,'_'),function (x) str_sub(x[3],2,4)),
#                 run= sapply(str_split(files_list,'_'),function (x) parse_number(x[4])),
#                 cutoff_prob= sapply(str_split(files_list,'_'),function (x) parse_number(x[5])),
#                 type=sapply(str_split(files_list,'_'),function (x) paste(x[6:8],collapse='_'))
#                 )
# files$PS <- sprintf('%0.2d', files$PS)
# files$type <- str_remove_all(files$type,'_NA')
# 
# files %>% filter(type=='modules.csv') %>% group_by(PS,scenario,cutoff_prob) %>% count(type) %>% print(n=Inf)

## @knitr END

# Deiversity regimes ------------------------------------------------------

## @knitr diversity_regime_description
if (!file.exists('/media/Data/PLOS_Biol/Results/regime_summary_data.csv')){
  experiments <- subset(all_experiments,PS!='18')
  regime_summary_data <- NULL
  for (i in 1:nrow(experiments)){
    ps <- experiments$PS[i]
    scenario <- experiments$scenario[i]
    exp <- experiments$exp[i]
    cutoff_prob <- experiments$cutoff_prob[i]
    run <- experiments$run[i]
    print(paste('PS: ',ps,' | Scenario: ',scenario,' | exp: ',exp, ' | run: ',run,sep=''))
    tmp <- get_data(parameter_space = ps, scenario, experiment = exp, run, cutoff_prob, use_sqlite = F, tables_to_get = 'summary_general')[[1]]
    regime_summary_data <- rbind(regime_summary_data,tmp)
  }
  write_csv(regime_summary_data, '/media/Data/PLOS_Biol/Results/regime_summary_data.csv')
} else {
  regime_summary_data <- read_csv('/media/Data/PLOS_Biol/Results/regime_summary_data.csv')
}

regime_summary_data %>%
  mutate(scenario=factor(scenario, levels=c('S','G','N'))) %>% 
  select(-year, -month, -n_infected) %>% 
  gather(variable, value, -pop_id, -time, -exp, -PS, -scenario, -run) %>% 
  group_by(pop_id, time, exp, PS, scenario, variable) %>%
  summarise(value_mean=mean(value), value_sd=sd(value)) %>% # Need to average across runs
  filter(variable %in% monitored_variables) %>%
  mutate(variable=factor(variable, levels=c('prevalence', 'meanMOI','n_alleles','n_circulating_genes','n_circulating_strains','n_total_bites'))) %>% 
  ggplot(aes(x=scenario, y=value_mean, fill=scenario))+
  geom_boxplot()+
  scale_fill_manual(values=scenario_cols)+
  facet_grid(variable~PS, scales='free',labeller = gg_labels)+
  labs(x='Scenario', y='Variable value')+
  mytheme_no_legend

regime_summary_data %>% 
  mutate(scenario=factor(scenario, levels=c('S','G','N'))) %>% 
  ggplot(aes(x=month, y=EIR, color=scenario))+
  geom_boxplot()+
  # stat_summary(aes(group=scenario), fun.y=mean, geom="point", size=3)+
  # stat_summary(aes(group=scenario), fun.y=mean, geom="line", size=0.5)+
  facet_wrap(~PS,scales='free', labeller = gg_labels)+
  scale_color_manual(values=scenario_cols)+
  labs(x='Month',y='EIR')+
  mytheme_no_legend

## @knitr END


# Edge weight distributions ------------------------------------------------------


# Fig. 2 ------------------------------------------------------------------
# Figure has 3 scenarios in high diversity, showing the following:
# Module examples, Relative persistence, Temporal Diversity, mFst.
# Number of repertoires per module is not so useful.

PS_for_figure <- '18'
module_results <- persistence_df <- temporal_diversity <- mFst <- c()
experiments <- subset(all_experiments, PS==PS_for_figure)
exp <- experiments$exp[1]
for (scenario in unique(experiments$scenario)){
  print(paste(PS_for_figure,scenario,exp, sep=' | '))
  
  f <- paste('/media/Data/PLOS_Biol/Results/PS',PS_for_figure,'_',scenario,'_E',exp,'_module_results.csv',sep='')
  x <- fread(f)
  module_results <- rbind(module_results,x)
  
  f <- paste('/media/Data/PLOS_Biol/Results/PS',PS_for_figure,'_',scenario,'_E',exp,'_persistence_df.csv',sep='')
  x <- fread(f)
  persistence_df <- rbind(persistence_df,x)
  
  x <- read_csv(paste('/media/Data/PLOS_Biol/Results/PS',PS_for_figure,'_',scenario,'_E',exp,'_temporal_diversity.csv',sep=''))
  temporal_diversity <- rbind(temporal_diversity,x)
  
  x <- read_csv(paste('/media/Data/PLOS_Biol/Results/PS',PS_for_figure,'_',scenario,'_E',exp,'_mFst.csv',sep=''))
  mFst <- rbind(mFst,x)
}

module_results$scenario <- factor(module_results$scenario, levels=c('S','G','N'))
module_results <- as.tibble(module_results) 
module_results %<>% mutate(PS=str_pad(PS, 2, 'left', '0'),exp=str_pad(exp, 3, 'left', '0'))
persistence_df$scenario <- factor(persistence_df$scenario, levels=c('S','G','N'))
persistence_df <- as.tibble(persistence_df) 
persistence_df %<>% mutate(PS=str_pad(PS, 2, 'left', '0'),exp=str_pad(exp, 3, 'left', '0'))
temporal_diversity$scenario <- factor(temporal_diversity$scenario, levels=c('S','G','N'))
mFst$scenario <- factor(mFst$scenario, levels=c('S','G','N'))


panel_A <- module_results %>% 
  filter(run==1) %>%
  filter(scenario=='S') %>% 
  distinct(module, layer) %>% 
  group_by(module) %>% summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module))+
  geom_rect(size=2, color=scenario_cols[1])+
  manuscript_theme+theme(axis.title = element_blank())
panel_B <- module_results %>% 
  filter(run==1) %>%
  filter(scenario=='G') %>% 
  distinct(module, layer) %>% 
  group_by(module) %>% summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module))+
  geom_rect(size=2, color=scenario_cols[2])+
  manuscript_theme+theme(axis.title = element_blank())
panel_C <- module_results %>% 
  filter(run==1) %>%
  filter(scenario=='N') %>% 
  distinct(module, layer) %>% 
  group_by(module) %>% summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module))+
  geom_rect(size=2, color=scenario_cols[3])+
  manuscript_theme+theme(axis.title = element_blank())
Fig <- plot_grid(panel_A,panel_B,panel_C, labels=c('A','B','C'), ncol=3, align='vh', label_size = 18, scale=0.95)
y.grob <- textGrob("Module ID", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
x.grob <- textGrob("Time (months)", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
Fig_2ABC <- grid.arrange(arrangeGrob(Fig, left = y.grob, bottom = x.grob))


panel_D <- persistence_df %>%
  ggplot()+
  geom_density(data=subset(persistence_df, type=='Repertoire' & scenario=='S'), aes(relative_persistence),fill='gray')+
  geom_density(data=subset(persistence_df, type=='Module' & scenario=='S'), aes(relative_persistence), fill=scenario_cols[1])+
  manuscript_theme+theme(axis.title = element_blank(),panel.grid = element_blank())
panel_E <- persistence_df %>%
  ggplot()+
  geom_density(data=subset(persistence_df, type=='Repertoire' & scenario=='G'), aes(relative_persistence),fill='gray')+
  geom_density(data=subset(persistence_df, type=='Module' & scenario=='G'), aes(relative_persistence), fill=scenario_cols[2])+
  manuscript_theme+theme(axis.title = element_blank(),panel.grid = element_blank())
panel_F <- persistence_df %>%
  ggplot()+
  geom_density(data=subset(persistence_df, type=='Repertoire' & scenario=='N'), aes(relative_persistence),fill='gray')+
  geom_density(data=subset(persistence_df, type=='Module' & scenario=='N'), aes(relative_persistence), fill=scenario_cols[3])+
  manuscript_theme+theme(axis.title = element_blank(),panel.grid = element_blank())

Fig <- plot_grid(panel_D,panel_E,panel_F, labels=c('D','E','F'), ncol=3, align='vh', label_size = 18, scale=0.95)
y.grob <- textGrob("Density", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
x.grob <- textGrob("Relative persistence", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
Fig_2DEF <- grid.arrange(arrangeGrob(Fig, left = y.grob, bottom = x.grob))

panel_D_inset <- persistence_df %>% 
  filter(scenario=='S') %>% 
  mutate(type=ifelse(type=='Module','M','R')) %>% 
  ggplot(aes(x=type, group=type, y=relative_persistence, fill=type))+
  geom_boxplot(outlier.size = 0, size=0.3)+
  facet_grid(~scenario)+
  scale_fill_manual(values = c(scenario_cols[1],'gray'))+
  labs(x='', y='Relative persistence')+
  manuscript_theme+theme(panel.grid = element_blank())
# panel_D <- panel_D+annotation_custom(grob=ggplotGrob(panel_D_inset), xmin=0.3,xmax=1,ymin=30,ymax=55) # High diversity
# panel_D <- panel_D+annotation_custom(grob=ggplotGrob(panel_D_inset), xmin=0.3,xmax=1,ymin=4,ymax=9) # Low diversity
# panel_D <- panel_D+annotation_custom(grob=ggplotGrob(panel_D_inset), xmin=0.3,xmax=1,ymin=25,ymax=50) # Medium diversity
panel_E_inset <- persistence_df %>% 
  filter(scenario=='G') %>% 
  mutate(type=ifelse(type=='Module','M','R')) %>% 
  ggplot(aes(x=type, group=type, y=relative_persistence, fill=type))+
  geom_boxplot(outlier.size = 0, size=0.3)+
  facet_grid(~scenario)+
  scale_fill_manual(values = c(scenario_cols[2],'gray'))+
  labs(x='', y='Relative persistence')+
  manuscript_theme+theme(panel.grid = element_blank())
# panel_E <- panel_E+annotation_custom(grob=ggplotGrob(panel_E_inset), xmin=0.3,xmax=1,ymin=30,ymax=55) # High diversity
# panel_E <- panel_E+annotation_custom(grob=ggplotGrob(panel_E_inset), xmin=0.3,xmax=1,ymin=4,ymax=9) # Low diversity
# panel_E <- panel_E+annotation_custom(grob=ggplotGrob(panel_E_inset), xmin=0.3,xmax=1,ymin=22,ymax=45) # Medium diversity
panel_F_inset <- persistence_df %>% 
  filter(scenario=='N') %>% 
  mutate(type=ifelse(type=='Module','M','R')) %>% 
  ggplot(aes(x=type, group=type, y=relative_persistence, fill=type))+
  geom_boxplot(outlier.size = 0, size=0.3)+
  facet_grid(~scenario)+
  scale_fill_manual(values = c(scenario_cols[3],'gray'))+
  labs(x='', y='Relative persistence')+
  manuscript_theme+theme(panel.grid = element_blank())
# panel_F <- panel_F+annotation_custom(grob=ggplotGrob(panel_F_inset), xmin=0.3,xmax=1,ymin=80,ymax=160) # High diversity
# panel_F <- panel_F+annotation_custom(grob=ggplotGrob(panel_F_inset), xmin=0.3,xmax=1,ymin=40,ymax=78) # Seasonality
# panel_F <- panel_F+annotation_custom(grob=ggplotGrob(panel_F_inset), xmin=0.3,xmax=1,ymin=7,ymax=16) # Low diversity
# panel_F <- panel_F+annotation_custom(grob=ggplotGrob(panel_F_inset), xmin=0.3,xmax=1,ymin=21,ymax=42) # Medium diversity
pdf('/home/shai/Dropbox/PLoS Biol/Fig_SI_seasonality_results_inset_F.pdf', 6,6)
svg('/home/shai/Dropbox/PLoS Biol/Fig_SI_seasonality_results_inset_D.svg', 6,6)
png('/home/shai/Dropbox/PLoS Biol/Fig_SI_seasonality_results_inset_F.png', 2000,2000, res=600)
panel_F_inset
dev.off()

dev.off()
png('/home/shai/Dropbox/PLoS Biol/Fig_structure_low.png', 4480*1.5,4490*1.5,units = 'px', res = 600)
pdf('/home/shai/Dropbox/PLoS Biol/Fig_structure_medium.pdf', 16,12)
pdf('/home/shai/Dropbox/PLoS Biol/Fig_SI_seasonality_results', 16,12)
# svg('/home/shai/Dropbox/PLoS Biol/Fig_structure_high.svg', 10.66667,8)
plot_grid(Fig_2ABC,Fig_2DEF, nrow=2, align='vh')
dev.off()

# pdf('/home/shai/Dropbox/Images for presentations/BGU Job talk/fig_2_panel_F.pdf',6,6)
# png('/home/shai/Dropbox/Images for presentations/BGU Job talk/fig_2_panel_A.png',1400,1400,res=200)
# panel_A
# dev.off()


# temporal_diversity$temp_div2 <- temporal_diversity$D/temporal_diversity$relative_persistence
# panel_G <- temporal_diversity %>%
#   filter(scenario=='S') %>%
#   ggplot()+
#   geom_density(aes(temp_div2), fill=scenario_cols[1])+
#   manuscript_theme+theme(axis.title = element_blank())
# panel_H <- temporal_diversity %>%
#   filter(scenario=='G') %>%
#   ggplot()+
#   geom_density(aes(temp_div2), fill=scenario_cols[2])+
#   manuscript_theme+theme(axis.title = element_blank())
# panel_I <- temporal_diversity %>%
#   filter(scenario=='N') %>%
#   ggplot()+
#   geom_density(aes(temp_div2), fill=scenario_cols[3])+
#   manuscript_theme+theme(axis.title = element_blank())
# Fig <- plot_grid(panel_G,panel_H,panel_I, labels=c('G','H','I'), ncol=3, align='vh', label_size = 18, scale=0.95)
# # y.grob <- textGrob("Density", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
# x.grob <- textGrob("Temporal diversity", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
# Fig_2GHI <- grid.arrange(arrangeGrob(Fig, left = y.grob, bottom = x.grob))
# 
# 
# panel_J <- mFst  %>% 
#   filter(scenario=='S') %>% 
#   ggplot()+
#   geom_density(aes(mFst), fill=scenario_cols[1])+
#   manuscript_theme+theme(axis.title = element_blank())
# panel_K <- mFst %>% 
#   filter(scenario=='G') %>% 
#   ggplot()+
#   geom_density(aes(mFst), fill=scenario_cols[2])+
#   manuscript_theme+theme(axis.title = element_blank())
# panel_L <- mFst %>% 
#   filter(scenario=='N') %>% 
#   ggplot()+
#   geom_density(aes(mFst), fill=scenario_cols[3])+
#   manuscript_theme+theme(axis.title = element_blank())
# Fig <- plot_grid(panel_J,panel_K,panel_L, labels=c('J','K','L'), ncol=3, align='vh', label_size = 18, scale=0.95)
# y.grob <- textGrob("Density", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
# x.grob <- textGrob("mFst", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
# Fig_2JKL <- grid.arrange(arrangeGrob(Fig, left = y.grob, bottom = x.grob))


# title <- ggdraw() + draw_label("High diversity", size=20)
# dev.off()
# dev.off()
# pdf('Results/Fig2_PS06.pdf', 16,12)
# plot_grid(title, Fig_2ABC,Fig_2DEF,Fig_2GHI, nrow=4, align='vh', rel_heights = c(0.096,0.32,0.32,0.32))
# # plot_grid(Fig_2JKL)
# dev.off()


# pdf('Results/Fig2_PS06_temporal_diversity.pdf', 16,10)
# temporal_diversity %>% 
#   ggplot(aes(x=scenario, group=scenario, y=statistic, fill=scenario))+
#   geom_boxplot(alpha=0.8)+
#   scale_fill_manual(values = scenario_cols)+
#   labs(x='', y='Temporal diversity', title='High diversity')+
#   manuscript_theme
# dev.off()

# Evenness in modules ----------------------------------------
# calculate_evenness_per_layer <- function(scen,r){
#   x <- module_results %>% 
#     filter(scenario==scen) %>% 
#     filter(run==r)
#   evenness <- NULL
#   for (l in unique(x$layer)){
#     x_l <- subset(x, layer==l)
#     d <- xtabs(~module+strain_id,data=x_l)
#     J <- vegan::diversity(d)/log(ncol(d))
#     tmp <- tibble(scenario=scen,run=r,layer=l, module=rownames(d), J=J)
#     evenness <- rbind(evenness, tmp)
#   }
#   return(evenness)
# }
# calculate_evenness <- function(scen,r){
#   x <- module_results %>% 
#     filter(scenario==scen) %>% 
#     filter(run==r)
#   d <- xtabs(~module+strain_id,data=x)
#   J <- vegan::diversity(d)/log(ncol(d))
#   return(tibble(scenario=scen,run=r,module=rownames(d),J=J))
# }



calculate_module_size_evenness_per_layer <- function(scen,r, module_results){
  x <- module_results %>% 
    filter(scenario==scen) %>% 
    filter(run==r)
  evenness <- NULL
  for (l in unique(x$layer)){
    x_l <- subset(x, layer==l)
    d <- x_l %>% group_by(module) %>% summarise(n=n())
    d <- d$n
    J <- vegan::diversity(d)/log(length(d))
    tmp <- tibble(scenario=scen,run=r,layer=l, n_modules=length(d), J=J)
    evenness <- rbind(evenness, tmp)
  }
  return(evenness)
}

module_size_evenness_per_layer <- NULL
for(i in 1:nrow(experiments)){
  scen <- experiments[i,'scenario']
  r <- experiments[i,'run']
  print(paste(scen,'|',r))
  y <- calculate_module_size_evenness_per_layer(scen,r,module_results = module_results)
  module_size_evenness_per_layer <- rbind(module_size_evenness_per_layer, y)
}

png('/home/shai/Dropbox/PLoS Biol/Fig_SI_evenness.png', 4480,4490, units = 'px', res = 600)
pdf('/home/shai/Dropbox/PLoS Biol/Fig_SI_evenness.pdf', 6,6)
module_size_evenness_per_layer %>% 
  mutate(scenario=factor(scenario, levels=c('S','G','N'))) %>% 
  ggplot(aes(x=J, y=n_modules, color=scenario))+
  geom_point(alpha=0.5)+
  scale_color_manual(values = scenario_cols)+
  labs(x='Evenness (J)', y='Number of modules')+
  mytheme_no_legend
dev.off()

# pdf('/home/shai/Dropbox/Images for presentations/BGU Job talk/Fig_SI_module_evenness.pdf',6,6)
# module_size_evenness_per_layer %>% 
#   mutate(scenario=factor(scenario, levels=c('S','G','N'))) %>% 
#   ggplot(aes(x=J, y=n_modules, color=scenario))+
#   geom_point(alpha=0.5)+
#   scale_color_manual(values = scenario_cols)+
#   labs(x='Evenness (J)', y='Number of modules')+
#   mytheme_no_legend
# dev.off()

pdf('/home/shai/Dropbox/PLoS Biol/Fig_SI_evenness_high.pdf', 6,6)
module_size_evenness_per_layer %>% 
  mutate(scenario=factor(scenario, levels=c('S','G','N'))) %>% 
  ggplot(aes(x=scenario, y=J, fill=scenario))+
  geom_boxplot()+
  scale_fill_manual(values = scenario_cols)+
  labs(x='Evenness (J)', y='Number of modules')+
  mytheme_no_legend
dev.off()



J_s <- module_size_evenness_per_layer %>% filter(scenario=='S')
J_g <- module_size_evenness_per_layer %>% filter(scenario=='G')
t.test(J_s$J,J_g$J)
wilcox.test(J_s$J,J_g$J)
median(J_s$J, na.rm = T)
median(J_g$J, na.rm = T)

# Seasonality --------------------------------------------
experiments <- subset(all_experiments,PS=='18')
regime_summary_data <- NULL
for (i in 1:nrow(experiments)){
  ps <- experiments$PS[i]
  scenario <- experiments$scenario[i]
  exp <- experiments$exp[i]
  cutoff_prob <- experiments$cutoff_prob[i]
  run <- experiments$run[i]
  print(paste('PS: ',ps,' | Scenario: ',scenario,' | exp: ',exp, ' | run: ',run,sep=''))
  tmp <- get_data(parameter_space = ps, scenario, experiment = exp, run, cutoff_prob, use_sqlite = F, tables_to_get = 'summary_general')[[1]]
  regime_summary_data <- rbind(regime_summary_data,tmp)
}

png('Results/Fig_SI_seasonality_variables.png', 1600,1200, res = 150)
regime_summary_data %>%
  mutate(scenario=factor(scenario, levels=c('S','G','N'))) %>% 
  select(-year, -n_infected) %>% 
  gather(variable, value, -pop_id, -time, -month, -exp, -PS, -scenario, -run) %>% 
  group_by(pop_id, month, time, exp, PS, scenario, variable) %>%
  summarise(value_mean=mean(value), value_sd=sd(value)) %>% # Need to average across runs
  filter(variable %in% monitored_variables) %>%
  mutate(variable=factor(variable, levels=c('prevalence', 'meanMOI','n_alleles','n_circulating_genes','n_circulating_strains','n_total_bites'))) %>% 
  ggplot(aes(x=month, y=value_mean, color=scenario))+
  geom_boxplot()+
  stat_summary(aes(group=scenario, color=scenario), fun.y=mean, geom="line", size=0.5)+
  scale_color_manual(values=scenario_cols)+
  facet_grid(variable~scenario, scales='free',labeller = gg_labels)+
  labs(x='Month', y='Variable value')+
  mytheme_no_legend+theme(axis.text.x = element_text(angle=-90))
dev.off()

png('/home/shai/Dropbox/PLoS Biol/SI_seasonality_EIR.png', 4480*1.5,4490*1.5,units = 'px', res = 600)
pdf('/home/shai/Dropbox/PLoS Biol/SI_seasonality_EIR.pdf', 16,12)
# png('/home/shai/Dropbox/Images for presentations/BGU Job talk/SI_seasonality_EIR.png',1400,1400,res=300)
regime_summary_data %>% 
  mutate(scenario=factor(scenario, levels=c('S','G','N'))) %>% 
  ggplot(aes(x=month, y=EIR, color=scenario))+
  geom_boxplot()+
  stat_summary(aes(group=scenario), fun.y=mean, geom="line", size=0.5)+
  scale_color_manual(values=scenario_cols)+
  labs(x='Month',y='EIR')+
  scale_y_continuous(breaks=seq(0,20,4))+
  facet_wrap(~scenario)+
  manuscript_theme
dev.off()




module_time_series <- module_results %>% 
  filter(scenario=='S') %>% 
  group_by(PS,scenario,run,layer) %>% 
  summarise(modules=length(unique(module)),
            repertoires=length(unique(strain_cluster))) %>% 
  group_by(PS,scenario,layer) %>% 
  summarise(modules_mean=mean(modules),
            repertoires_mean=mean(repertoires)) 
module_ts_modules <- ts(module_time_series$modules_mean, start=c(1,1), deltat = 1/12)
module_ts_repertoires <- ts(module_time_series$repertoires_mean, start=c(1,1), deltat = 1/12)
decomposedRes_modules <- stats::decompose(module_ts_modules, type="mult") # use type = "additive" for additive components
decomposedRes_repertoires <- stats::decompose(module_ts_repertoires, type="mult") # use type = "additive" for additive components

dev.off()
# pdf('Results/Fig2_PS18_time_series.pdf', 16,12)
png('Results/Fig2_PS18_time_series_repertoires.png', 1600,1200, res=150)
module_results %>% 
  group_by(PS,scenario,run,layer) %>% 
  summarise(repertoires=length(unique(strain_cluster))) %>% 
  group_by(PS,scenario,layer) %>% 
  summarise(repertoires_mean=mean(repertoires),
            repertoires_sd=sd(repertoires)) %>% 
  ggplot()+
  geom_errorbar(aes(x=layer, ymax=repertoires_mean+repertoires_sd, ymin=repertoires_mean-repertoires_sd, group=scenario),color='gray50')+
  geom_line(aes(x=layer, y=repertoires_mean, color=scenario))+
  scale_color_manual(values = scenario_cols)+
  labs(y= '# repertoires', x='Time (months)', title='Repertoires time series')+mytheme
dev.off()

png('Results/Fig2_PS18_time_series_modules.png', 1600,1200, res=150)
module_results %>% 
  group_by(PS,scenario,run,layer) %>% 
  summarise(modules=length(unique(module))) %>% 
  group_by(PS,scenario,layer) %>% 
  summarise(modules_mean=mean(modules),
            modules_sd=sd(modules)) %>% 
  ggplot()+
  geom_errorbar(aes(x=layer, ymax=modules_mean+modules_sd, ymin=modules_mean-modules_sd, group=scenario),color='gray50')+
  geom_line(aes(x=layer, y=modules_mean, color=scenario))+
  scale_color_manual(values = scenario_cols)+
  labs(y= '# modules', x='Time (months)', title='Modules time series')+mytheme
dev.off()

ccf(module_ts_modules,module_ts_repertoires, main='CCF between # repertoires and # modules')
plot (decomposedRes_repertoires) # see plot below
plot(decomposedRes_modules) # see plot below
# dev.off()


# Cutoff plots (Selection) ------------------------------------------------

results_cutoff_S <- fread('/media/Data/PLOS_Biol/Results/results_cutoff_S.csv', 
                        colClasses = list(integer=c(1,2,3,4,6,10),
                                          character=c(5,7,8,9),
                                          double=11))
results_cutoff_S <- as.tibble(results_cutoff_S)

experiments <- results_cutoff_S %>% distinct(scenario,PS,cutoff_prob,run)
scen <- 'S'
module_size_evenness_per_layer_cutoff <- NULL
for(i in 1:nrow(experiments)){
  ps <- experiments$PS[i]
  r <- experiments$run[i]
  cutoff <- experiments$cutoff_prob[i]
  print(paste(i,ps,cutoff,r,sep=' | '))
  x <- results_cutoff_S %>% 
    filter(PS==ps, cutoff_prob==cutoff) %>% 
    filter(run==r)
  evenness <- NULL
  for (l in unique(x$layer)){
    x_l <- subset(x, layer==l)
    d <- x_l %>% group_by(module) %>% summarise(n=n())
    d <- d$n
    J <- vegan::diversity(d)/log(length(d))
    tmp <- tibble(scenario=scen,run=r,layer=l, n_modules=length(d), J=J)
    evenness <- rbind(evenness, tmp)
  }
  evenness$PS <- ps
  evenness$run <- run
  evenness$cutoff_prob <- cutoff
  module_size_evenness_per_layer_cutoff <- rbind(module_size_evenness_per_layer_cutoff, evenness)
}

module_size_evenness_per_layer_cutoff %>% 
  ggplot(aes(group=cutoff_prob, x=cutoff_prob, y=J, fill=PS))+
  geom_boxplot()+
  scale_fill_manual(values = ps_cols)+
  labs(y='Evenness (J)', x='Cutoff')+
  facet_wrap(~PS)+
  mytheme

persistence_df_cutoff_S <- read_csv('/media/Data/PLOS_Biol/Results/persistence_df_cutoff_S.csv', col_types = 'cidciiidc')
temporal_diversity_cutoff_S <- read_csv('/media/Data/PLOS_Biol/Results/temporal_diversity_cutoff_S.csv', col_types = 'ccidiiiiddd')

## Relative persistence
panel_A <- persistence_df_cutoff_S %>% 
  filter(PS=='04', type=='Module') %>% 
  filter(cutoff_prob>=0.25) %>% 
  ggplot(aes(x=cutoff_prob, y=relative_persistence, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[1])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
  scale_x_continuous(breaks = seq(0.05,0.95,0.15))+
  manuscript_theme+theme(axis.title = element_blank())
panel_B <- persistence_df_cutoff_S %>% 
  filter(PS=='05', type=='Module') %>% 
  filter(cutoff_prob>=0.25) %>% 
  ggplot(aes(x=cutoff_prob, y=relative_persistence, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[2])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
  scale_x_continuous(breaks = seq(0.05,0.95,0.15))+
  manuscript_theme+theme(axis.title = element_blank())
panel_C <- persistence_df_cutoff_S %>% 
  filter(PS=='06', type=='Module') %>% 
  filter(cutoff_prob>=0.25) %>% 
  ggplot(aes(x=cutoff_prob, y=relative_persistence, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[3])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
  scale_x_continuous(breaks = seq(0.05,0.95,0.15))+
  manuscript_theme+theme(axis.title = element_blank())
Fig <- plot_grid(panel_A,panel_B,panel_C, labels=c('A','B','C'), ncol=3, align='vh', label_size = 18, scale=0.95)
y.grob <- textGrob("Relative persistence", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
# x.grob <- textGrob("Cut off", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
Fig_cutoff_ABC <- grid.arrange(arrangeGrob(Fig, left = y.grob))

## Modules per layer
panel_D <- results_cutoff_S %>% 
  filter(PS=='04') %>% 
  filter(cutoff_prob>=0.25) %>% 
  group_by(scenario, run, cutoff_prob, layer) %>% 
  summarise(modules_per_layer=length(unique(module))) %>% 
  ggplot(aes(x=cutoff_prob, y=modules_per_layer, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[1])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.15))+
  manuscript_theme+theme(axis.title = element_blank())
panel_E <- results_cutoff_S %>% 
  filter(PS=='05') %>% 
  filter(cutoff_prob>=0.25) %>% 
  group_by(scenario, run, cutoff_prob, layer) %>% 
  summarise(modules_per_layer=length(unique(module))) %>% 
  ggplot(aes(x=cutoff_prob, y=modules_per_layer, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[2])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.15))+
  manuscript_theme+theme(axis.title = element_blank())
panel_F <- results_cutoff_S %>% 
  filter(PS=='06') %>% 
  filter(cutoff_prob>=0.25) %>% 
  group_by(scenario, run, cutoff_prob, layer) %>% 
  summarise(modules_per_layer=length(unique(module))) %>% 
  ggplot(aes(x=cutoff_prob, y=modules_per_layer, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[3])+
  stat_summary(fun.y=mean, color="red", geom="point", shape=18, size=3,show.legend = FALSE)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.15))+
  manuscript_theme+theme(axis.title = element_blank())
Fig <- plot_grid(panel_D,panel_E,panel_F, labels=c('D','E','F'), ncol=3, align='vh', label_size = 18, scale=0.95)
y.grob <- textGrob("Modules per layer", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
# x.grob <- textGrob("Cut off", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
Fig_cutoff_DEF <- grid.arrange(arrangeGrob(Fig, left = y.grob))


## Evenness in module size
panel_G <- module_size_evenness_per_layer_cutoff %>%
  filter(PS=='04') %>% 
  ggplot(aes(x=cutoff_prob, y=J, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[1])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
  scale_x_continuous(breaks = seq(0.05,0.95,0.15))+
  manuscript_theme+theme(axis.title = element_blank())
panel_H <- module_size_evenness_per_layer_cutoff %>%
  filter(PS=='05') %>% 
  ggplot(aes(x=cutoff_prob, y=J, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[2])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
  scale_x_continuous(breaks = seq(0.05,0.95,0.15))+
  manuscript_theme+theme(axis.title = element_blank())
panel_I <- module_size_evenness_per_layer_cutoff %>%
  filter(PS=='06') %>% 
  ggplot(aes(x=cutoff_prob, y=J, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[3])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
  scale_x_continuous(breaks = seq(0.05,0.95,0.15))+
  manuscript_theme+theme(axis.title = element_blank())
Fig <- plot_grid(panel_G,panel_H,panel_I, labels=c('G','H','I'), ncol=3, align='vh', label_size = 18, scale=0.95)
y.grob <- textGrob("Evenness (J)", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
x.grob <- textGrob("Cut off", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
Fig_cutoff_GHI <- grid.arrange(arrangeGrob(Fig, left = y.grob, bottom = x.grob))

png('/home/shai/Dropbox/PLoS Biol/Fig_SI_cutoff_sweep_edge_dist.png', 4480*1.5,4490*1.5,units = 'px', res = 600)
pdf('/home/shai/Dropbox/PLoS Biol/Fig_SI_cutoff_sweep_edge_dist.pdf', 16,12)
# svg('/home/shai/Dropbox/PLoS Biol/Fig_structure_high.svg', 10.66667,8)
plot_grid(Fig_cutoff_ABC,Fig_cutoff_DEF,Fig_cutoff_GHI, nrow=3, align='vh')
dev.off()

## Temporal diveristy
# panel_G <- temporal_diversity_cutoff_S %>%
#   filter(PS=='04') %>% 
#   ggplot(aes(x=cutoff_prob, y=statistic, group=cutoff_prob))+
#   geom_boxplot(outlier.size = 0, fill=ps_cols[1])+
#   stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
#   scale_x_continuous(breaks = seq(0.05,0.95,0.15))+
#   manuscript_theme+theme(axis.title = element_blank())
# panel_H <- temporal_diversity_cutoff_S %>%
#   filter(PS=='05') %>% 
#   ggplot(aes(x=cutoff_prob, y=statistic, group=cutoff_prob))+
#   geom_boxplot(outlier.size = 0, fill=ps_cols[2])+
#   stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
#   scale_x_continuous(breaks = seq(0.05,0.95,0.15))+
#   manuscript_theme+theme(axis.title = element_blank())
# panel_I <- temporal_diversity_cutoff_S %>%
#   filter(PS=='06') %>% 
#   ggplot(aes(x=cutoff_prob, y=statistic, group=cutoff_prob))+
#   geom_boxplot(outlier.size = 0, fill=ps_cols[3])+
#   stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
#   scale_x_continuous(breaks = seq(0.05,0.95,0.15))+
#   manuscript_theme+theme(axis.title = element_blank())
# Fig <- plot_grid(panel_G,panel_H,panel_I, labels=c('G','H','I'), ncol=3, align='vh', label_size = 18, scale=0.95)
# y.grob <- textGrob("Temporal diversity", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
# x.grob <- textGrob("Cut off", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
# Fig_cutoff_GHI <- grid.arrange(arrangeGrob(Fig, left = y.grob, bottom = x.grob))

# title <- ggdraw() + draw_label("Cutoff results", size=20)
# plot_grid(title, Fig_cutoff_ABC,Fig_cutoff_DEF,Fig_cutoff_GHI, nrow=4, align='vh', rel_heights = c(0.096,0.32,0.32,0.32))
pdf('Results/Fig_SI_cutoff_S.pdf', 16,12)
plot_grid(title, Fig_cutoff_ABC,Fig_cutoff_DEF,Fig_cutoff_GHI, nrow=4, align='vh', rel_heights = c(0.096,0.32,0.32,0.32))
dev.off()


# Also plot an example of modules
results_cutoff_S %>% 
  # filter(PS=='06') %>% 
  filter(run==1) %>%
  filter(cutoff_prob>=0.5) %>% 
  distinct(PS,cutoff_prob, module, layer) %>% 
  group_by(PS,cutoff_prob,module) %>% summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module, color=PS, fill=PS))+
  geom_rect(size=1)+
  scale_color_manual(values=ps_cols)+
  scale_fill_manual(values=ps_cols)+
  facet_grid(cutoff_prob~PS, scales = 'free')+
  mytheme+theme(axis.title = element_blank(), legend.position = 'none')
dev.off()


png('Results/Fig_SI_cutoff_S.png', 1600,1200,res=150)
plot_grid(Fig_cutoff_ABC,Fig_cutoff_DEF,Fig_cutoff_GHI, nrow=3, align='vh')
dev.off()
png('Results/Fig_SI_cutoff_S_modules.png', 1600,1200,res=150)
results_cutoff_S %>% 
  # filter(PS=='06') %>% 
  filter(run==1) %>%
  filter(cutoff_prob>=0.5) %>% 
  distinct(PS,cutoff_prob, module, layer) %>% 
  group_by(PS,cutoff_prob,module) %>% summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module, color=PS, fill=PS))+
  geom_rect(size=1)+
  scale_color_manual(values=ps_cols)+
  scale_fill_manual(values=ps_cols)+
  facet_grid(cutoff_prob~PS, scales = 'free')+
  mytheme+theme(axis.title = element_blank(), legend.position = 'none')
dev.off()




# Cutoff plots (all scenarios) --------------------------------------------
results_cutoff_N <- fread('/media/Data/PLOS_Biol/Results/results_cutoff_N.csv', 
                        colClasses = list(integer=c(1,2,3,4,6,10),
                                          character=c(5,7,8,9),
                                          double=11))
results_cutoff_N <- as.tibble(results_cutoff_N)
persistence_df_cutoff_N <- read_csv('/media/Data/PLOS_Biol/Results/persistence_df_cutoff_N.csv', col_types = 'cidciiidc')
temporal_diversity_cutoff_N <- read_csv('/media/Data/PLOS_Biol/Results/temporal_diversity_cutoff_N.csv', col_types = 'ccidiiiiddd')

results_cutoff_G <- fread('/media/Data/PLOS_Biol/Results/results_cutoff_G.csv', 
                          colClasses = list(integer=c(1,2,3,4,6,10),
                                            character=c(5,7,8,9),
                                            double=11))
results_cutoff_G <- as.tibble(results_cutoff_G)
persistence_df_cutoff_G <- read_csv('/media/Data/PLOS_Biol/Results/persistence_df_cutoff_G.csv', col_types = 'cidciiidc')
temporal_diversity_cutoff_G <- read_csv('/media/Data/PLOS_Biol/Results/temporal_diversity_cutoff_G.csv', col_types = 'ccidiiiiddd')

## Join the scenarios to one data frame
results_cutoff_SGN <- rbind(results_cutoff_S,results_cutoff_G,results_cutoff_N)
results_cutoff_SGN$scenario <- factor(results_cutoff_SGN$scenario, levels=c('S','G','N'))

persistence_df_cutoff_S$scenario <- 'S'
persistence_df_cutoff_G$scenario <- 'G'
persistence_df_cutoff_N$scenario <- 'N'
persistence_df_cutoff_SGN <- rbind(persistence_df_cutoff_S,persistence_df_cutoff_G,persistence_df_cutoff_N)
persistence_df_cutoff_SGN$scenario <- factor(persistence_df_cutoff_SGN$scenario, levels=c('S','G','N'))

temporal_diversity_cutoff_SGN <- rbind(temporal_diversity_cutoff_S,temporal_diversity_cutoff_G,temporal_diversity_cutoff_N)
temporal_diversity_cutoff_SGN$scenario <- factor(temporal_diversity_cutoff_SGN$scenario, levels=c('S','G','N'))

## Plotting

pdf('Results/Fig_SI_cutoff_SNG.pdf', 16,12)

persistence_df_cutoff_SGN %>% 
  filter(type=='Module') %>% 
  group_by(scenario, PS, cutoff_prob) %>% 
  summarise(mean_relative_persistence=mean(relative_persistence),
            sd_relative_persistence=sd(relative_persistence),
            median_relative_persistence=median(relative_persistence)
  ) %>% ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_relative_persistence))+geom_line(aes(y=mean_relative_persistence))+
  # geom_errorbar(aes(x=cutoff_prob, ymin=mean_relative_persistence-sd_relative_persistence, ymax=mean_relative_persistence+sd_relative_persistence),color='gray')+
  geom_point(aes(y=median_relative_persistence), shape=15)+geom_line(aes(y=median_relative_persistence), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median module persistence')+
  scale_color_manual(values=scenario_cols)+mytheme


## Modules per layer
results_cutoff_SGN %>% 
  group_by(PS, scenario, run, cutoff_prob, layer) %>% 
  summarise(modules_per_layer=length(unique(module))) %>% 
  group_by(PS, scenario, cutoff_prob) %>%
  summarise(mean_modules_per_layer=mean(modules_per_layer),
            median_modules_per_layer=median(modules_per_layer)) %>%
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_modules_per_layer))+geom_line(aes(y=mean_modules_per_layer))+
  geom_point(aes(y=median_modules_per_layer), shape=15)+geom_line(aes(y=median_modules_per_layer), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median modules per layer')+
  scale_color_manual(values=scenario_cols)+mytheme

# Temporal diversity
x <- temporal_diversity_cutoff_SGN %>% 
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

x %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_statistic))+geom_line(aes(y=mean_statistic))+
  geom_point(aes(y=median_statistic), shape=15)+geom_line(aes(y=median_statistic), linetype='dashed')+
  # geom_errorbar(aes(ymin=ymin_statistic_ci,ymax=ymax_statistic_ci),alpha=0.6, width=0.01, size=1)+
  facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median Temporal Diversity')+
  scale_color_manual(values=scenario_cols)+mytheme


dev.off()



# Sensitivity analysis ----------------------------------------------------

sensitivity_experiments <- expand.grid(PS=as.character(100:211),
                               scenario=c('S','G'), 
                               run=1,
                               cutoff_prob=0.85,
                               stringsAsFactors = F)

run <- 1
cutoff_prob <- 0.85
module_results_sensitivity <- c()
temporal_diversity_sensitivity <- NULL
exp='001'
for (i in 1:nrow(sensitivity_experiments)){
  PS <- sensitivity_experiments$PS[i]
  scenario <- sensitivity_experiments$scenario[i]
  x <- get_modularity_results(PS,scenario,exp,run,cutoff_prob, folder = '/media/Data/PLOS_Biol/Results/sensitivity_analysis/')
  module_results_sensitivity <- rbind(module_results_sensitivity,x)
  x <- get_temporal_diversity(PS,scenario,exp,run,cutoff_prob, folder = '/media/Data/PLOS_Biol/Results/sensitivity_analysis/')
  temporal_diversity_sensitivity <- rbind(temporal_diversity_sensitivity, x)
}
module_results_sensitivity$strain_cluster <- as.integer(module_results_sensitivity$strain_cluster)
module_results_sensitivity <- as.tibble(module_results_sensitivity)
temporal_diversity_sensitivity <- as.tibble(temporal_diversity_sensitivity)
# write_csv(module_results_sensitivity, '/media/Data/PLOS_Biol/Results/module_results_sensitivity.csv')

x <- module_results_sensitivity %>%
  filter(layer<=300) %>% 
  mutate(scenario=factor(scenario, levels=c('S','G'))) %>% 
  select(scenario, PS, run, layer, module) %>% 
  group_by(scenario, PS,run,module) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) 

panel_A <- x %>% 
  ggplot(aes(x=relative_persistence, fill=scenario))+
  geom_density()+
  scale_fill_manual(values = c('red','blue'))+
  labs(x='Relativer persistence', y='Density')+
  manuscript_theme+theme(axis.text = element_text(size=20),
                         axis.title = element_text(size=20))
panel_B <- temporal_diversity_sensitivity %>% 
  mutate(scenario=factor(scenario, levels=c('S','G'))) %>% 
  ggplot(aes(x=scenario, y=D_normalized, fill=scenario))+
  geom_boxplot()+
  scale_fill_manual(values = c('red','blue'))+
  labs(x='Scenario', y='Evenness (J)')+
  manuscript_theme+theme(axis.text = element_text(size=20),
                         axis.title = element_text(size=20))
Fig <- plot_grid(panel_A,panel_B, labels=c('A','B'), ncol=2, align='vh', label_size = 18)
pdf('~/Dropbox/PLoS Biol/Fig_SI_sensitivity.pdf',16,9)
Fig
dev.off()
# # Module examples Can only be plotted for one specific run There may be gaps in
# # the layers for a given module because modules are based on sampled
# # layer-repertoire tuples, and repertoires may not be sampled in every layer in the ABM.
# # So there are two options to plot: one that shows the gaps, and one that does not.
# module_results_sensitivity %>% 
#   filter(PS==150) %>%
#   distinct(module, layer, PS, scenario) %>% 
#   group_by(scenario,module) %>% summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
#   ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module, color=scenario))+
#   geom_rect(size=2)+
#   scale_color_manual(values = scenario_cols)+
#   labs(y= 'module ID', x='Time (months)')+
#   facet_wrap(~scenario, scales='free', labeller = gg_labels)+manuscript_theme

# Relative persistence
module_persistence_sensitivity <- module_results_sensitivity %>% 
  select(scenario, PS, run, cutoff_prob, layer, module) %>% 
  group_by(scenario, PS,run,cutoff_prob,module) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Module') %>% 
  rename(id=module) %>% mutate(id=as.character(id))

strain_persistence_sensitivity <- module_results_sensitivity %>% 
  select(scenario, PS, run, cutoff_prob, layer, strain_cluster) %>% 
  group_by(scenario, PS,run,cutoff_prob,strain_cluster) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Repertoire') %>% 
  rename(id=strain_cluster) %>% mutate(id=as.character(id))

persistence_df_sensitivity <- module_persistence_sensitivity %>% bind_rows(strain_persistence_sensitivity)
write_csv(persistence_df_sensitivity, '/media/Data/PLOS_Biol/Results/persistence_df_sensitivity.csv')


# Read sensitivity results and plot
module_results_sensitivity <- read_csv('/media/Data/PLOS_Biol/Results/module_results_sensitivity.csv')
module_results_sensitivity$PS <- as.character(module_results_sensitivity$PS)
persistence_df_sensitivity <- read_csv('/media/Data/PLOS_Biol/Results/persistence_df_sensitivity.csv')
persistence_df_sensitivity$PS <- as.character(persistence_df_sensitivity$PS )

persistence_df <- fread('/media/Data/PLOS_Biol/Results/PS06_S_E001_persistence_df.csv') # get the PS06 persistence results
persistence_df <- as.tibble(persistence_df)
persistence_df$PS <- str_pad(persistence_df$PS, width = 2, side = 'left', pad = '0')
persistence_df_sensitivity$id <- as.integer(persistence_df_sensitivity$id)
persistence_df$id <- as.integer(persistence_df$id)


persistence_df_sensitivity$PS <- as.character(persistence_df_sensitivity$PS)
persistence_df_sensitivity$id <- as.integer(persistence_df_sensitivity$id)


persistence_df_sensitivity %<>% 
  bind_rows(persistence_df) %>%
  filter(scenario=='S') %>%
  mutate(grp=ifelse(PS=='06','Main result','Sensitivity'))
# Compare module persistence of sensitivity to observed

pdf('~/Dropbox/PLoS Biol/Fig_SI_sensitivity.pdf',12,12)
persistence_df_sensitivity %>% 
  filter(type=='Module') %>% 
  # filter(PS %in% c('06',as.character(172))) %>% 
  ggplot()+
  geom_boxplot(aes(x=grp, y=relative_persistence, fill=grp), notch = F)+
  # geom_density(aes(x=relative_persistence, color=g))+
  # geom_density(data=subset(persistence_df_sensitivity, PS =='06'), aes(relative_persistence), fill=scenario_cols[1], alpha=0.5)+
  labs(x='', y='Relative persistence')+
  manuscript_theme
dev.off()

persistence_df_sensitivity %>% 
  filter(type=='Module') %>% 
  group_by(grp,PS,run) %>% 
  summarise(m=mean(relative_persistence)) %>% 
  # filter(PS %in% c('06',as.character(172))) %>% 
  ggplot()+
  geom_boxplot(aes(x=grp, y=m, fill=grp), notch = F)+
  # geom_density(aes(x=relative_persistence, color=g))+
  # geom_density(data=subset(persistence_df_sensitivity, PS =='06'), aes(relative_persistence), fill=scenario_cols[1], alpha=0.5)+
  labs(x='', y='Relative persistence')+
  manuscript_theme



names(temporal_diversity)[10] <- 'J'
names(temporal_diversity_sensitivity)[11] <- 'J'
temporal_diversity_sensitivity$PS <- as.character(temporal_diversity_sensitivity$PS)
temporal_diversity$grp <- 'Main' 
temporal_diversity_sensitivity$grp <- 'Sensitivity' 

x <- temporal_diversity_sensitivity %>% select(PS,module,relative_persistence,J,grp) %>% 
  bind_rows(select(.data = temporal_diversity, PS,module,relative_persistence,J,grp))

unique(x$grp)

x %>% 
  # ggplot(aes(x=grp,y=J,fill=grp))+geom_boxplot()
  ggplot(aes(x=J,fill=grp))+geom_density()


x <- module_results %>% 
  group_by(layer,run) %>% 
  summarise(n_modules=length(unique(module))) %>% 
  group_by(layer) %>% summarise(n_modules_mean=mean(n_modules))
y <- module_results_sensitivity %>% 
  # filter(PS=='211') %>% 
  group_by(layer,PS) %>% 
  summarise(n_modules=length(unique(module))) %>% 
  group_by(layer) %>% summarise(n_modules_mean=mean(n_modules))
x$grp <- 'M'
y$grp <- 'Sen'
x %>% bind_rows(y) %>% 
  ggplot(aes(x=layer, y=n_modules_mean, color=grp))+geom_line()
  ggplot(aes(x=layer, y=n_modules_mean))+geom_line()
# persistence_df_sensitivity %>% 
#   filter(type=='Repertoire') %>% 
#   ggplot()+
#   geom_density(aes(x=relative_persistence, fill=grp))+
#   labs(x='', y='Relative persistence')+
#   mytheme

mean(subset(persistence_df_sensitivity, grp=='Sensitivity' & type=='Module')$relative_persistence)
mean(subset(persistence_df_sensitivity, grp=='Main result' & type=='Module')$relative_persistence)
median(subset(persistence_df_sensitivity, grp=='Sensitivity' & type=='Module')$relative_persistence)
median(subset(persistence_df_sensitivity, grp=='Main result' & type=='Module')$relative_persistence)

wilcox.test(x = subset(persistence_df_sensitivity, grp=='Sensitivity' & type=='Repertoire')$relative_persistence,
            y = subset(persistence_df_sensitivity, grp=='Main result' & type=='Repertoire')$relative_persistence)

# Temporal diversity
temporal_diversity_sensitivity <- c()
scenario <- 'S'
run <- 1
cutoff_prob <- 0.85
for (i in 1:nrow(sensitivity_experiments)){
  PS <- sensitivity_experiments$PS[i]
  x <- get_temporal_diversity(PS = PS,scenario = scenario,exp = '001', run = run,cutoff_prob = cutoff_prob,folder = '/media/Data/PLOS_Biol/Results/sensitivity_analysis/')
  temporal_diversity_sensitivity <- rbind(temporal_diversity_sensitivity,x)
}
temporal_diversity_sensitivity$PS <- as.character(temporal_diversity_sensitivity$PS )
temporal_diversity_sensitivity <- as.tibble(temporal_diversity_sensitivity)

temporal_diversity <- fread('/media/Data/PLOS_Biol/Results/PS06_S_E001_temporal_diversity.csv') # get the PS06 persistence results
temporal_diversity <- as.tibble(temporal_diversity)
temporal_diversity$PS <- str_pad(temporal_diversity$PS, width = 2, side = 'left', pad = '0')

names(temporal_diversity)[10] <- 'D_normalized'
temporal_diversity_sensitivity %<>% 
  bind_rows(temporal_diversity) %>% 
  filter(scenario=='S') %>% 
  mutate(grp=ifelse(PS=='06','Main result','Sensitivity'))
panel_B <- temporal_diversity_sensitivity %>% 
  ggplot()+
  geom_boxplot(aes(x=grp, y=D_normalized, fill=grp), notch = T)+
  # geom_density(data=subset(temporal_diversity_sensitivity, PS !='06'), aes(statistic), fill='#1E9B95', alpha=0.5)+
  # geom_density(data=subset(temporal_diversity_sensitivity, PS =='06'), aes(statistic), fill=scenario_cols[1], alpha=0.5)+
  labs(x='', y='Temporal diversity')+
  manuscript_theme
# panel_B <- panel_B + draw_label('Main result in red; \nsensitivity simulations in teal', 0.2, 7, hjust = 0, vjust = 0, size = 20)

pdf('Results/Fig_SI_sesitivity_selection.pdf', 16,12)
plot_grid(panel_A,panel_B, labels = c('A','B'))
dev.off()

png('Results/Fig_SI_sesitivity_selection.png', 1600,1200, res=150)
plot_grid(panel_A,panel_B, labels = c('A','B'))
dev.off()


# 
# # mFst
# mFst_sensitivity <- c()
# scenario <- 'S'
# run <- 1
# cutoff_prob <- 0.85
# for (i in 1:nrow(sensitivity_experiments)){
#   PS <- sensitivity_experiments$PS[i]
#   x <- get_mFst(PS,scenario,run,cutoff_prob,folder = '/media/Data/PLOS_Biol/Results/sensitivity_analysis/')
#   mFst_sensitivity <- rbind(mFst_sensitivity,x)
# }
# mFst_sensitivity <- as.tibble(mFst_sensitivity)
# x <- mFst_sensitivity %>% 
#   bind_rows(mFst) %>% 
#   filter(scenario=='S') %>% 
#   mutate(grp=ifelse(PS=='06','Obs','Sens'))
# x %>% 
#   ggplot()+
#   geom_density(data=subset(x, PS !='06'), aes(mFst), fill='#1E9B95', alpha=0.4)+
#   geom_density(data=subset(x, PS =='06'), aes(mFst), fill=scenario_cols[1], alpha=0.4)+
#   labs(x='mFst', y='Density')+
#   mytheme


# Statistical analysis for differences in persistence ----------
persistence_df %>% 
  group_by(scenario, type) %>% 
  summarise(mean_relative_persistence=mean(relative_persistence),
            sd_relative_persistence=sd(relative_persistence),
            median_relative_persistence=median(relative_persistence))

for (scen in c('S','G','N')){
  print(scen)
  module_pers <- subset(persistence_df, type=='Module' & scenario==scen)
  rep_pers <- subset(persistence_df, type=='Repertoire' & scenario==scen)
  print(wilcox.test(module_pers$relative_persistence, rep_pers$relative_persistence))
}

# Repertoires per module
# module_results %>% 
#   mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% 
#   group_by(scenario, run, module) %>% 
#   summarise(repertoires_per_module=length(unique(strain_cluster))) %>% 
#   select(scenario, run, repertoires_per_module) %>% 
#   ggplot(aes(x=scenario, y=repertoires_per_module, fill=scenario))+
#   geom_boxplot()+
#   scale_fill_manual(values=scenario_cols)+mytheme
# 
# x=module_results %>%
#   filter(scenario=='S') %>%
#   group_by(run, module) %>%
#   summarise(repertoires_per_module=length(unique(strain_cluster))) %>%
#   select(run,repertoires_per_module)
# y=module_results %>%
#   filter(scenario=='N') %>%
#   group_by(run, module) %>%
#   summarise(repertoires_per_module=length(unique(strain_cluster))) %>%
#   select(run,repertoires_per_module)
# t.test(x$repertoires_per_module,y$repertoires_per_module)
# wilcox.test(x$repertoires_per_module,y$repertoires_per_module)


# Epidemiological simulations ---------------------------------------------

get_epidemiology_results <- function(PS,scenario,exp,run,cutoff_prob,folder='/media/Data/PLOS_Biol/Results/epi/'){
  files <- list.files(path = folder, pattern = paste('PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,sep=''), full.names = T)
  if (length(files)==0){
    print(paste('No files found for ','PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,sep=''))
    return(NULL)
  }
  res_df <- NULL
  for (f in files){
    print(f)
    x <- read.csv(f)
    head(x)
  }
  res_df <- rbind(res_df, x)
  return(as.tibble(res_df))
}

experiments <- expand.grid(PS=c('06'),
                                scen=c('S','G','N'),
                                exp='001',
                                run=1:50,
                                stringsAsFactors = F)
cutoff_df <- tibble(PS=c('04','05','06','18'),
                    cutoff_prob=c(0.3,0.6,0.85,0.85))
experiments %<>% left_join(cutoff_df)

epi_results <- NULL
for (i in 1:nrow(experiments)){
  PS <- experiments$PS[i]
  scenario <- experiments$scen[i]
  cutoff_prob <- experiments$cutoff_prob[i]
  run <- experiments$run[i]
  tmp <- get_epidemiology_results(PS = PS,scenario = scenario,exp = '001',run = run, cutoff_prob = cutoff_prob, 
                                  folder='/media/Data/PLOS_Biol/Results/epi_alleles/')
  epi_results <- rbind(epi_results, tmp)
}

epi_results$PS <- str_pad(epi_results$PS, 2, 'left', '0')
epi_results$exp <- str_pad(epi_results$exp, 3, 'left', '0')

epi_results %>% group_by(scenario,PS,case) %>% summarise(n=length(unique(run))) %>% print(n=Inf)

# Limit the infections
infection_limit <- epi_results %>%
  group_by(PS,scenario,case) %>% 
  summarise(max_infections = max(infection_id)) %>% 
  group_by(PS,scenario) %>% 
  summarise(infection_limit=min(max_infections,na.rm = T))

make_panel_fig4 <- function(ps,scen,col){
  panel <- epi_results %>% 
    filter(PS==ps) %>% 
    filter(scenario==scen) %>% 
    group_by(case,infection_id) %>% 
    summarise(mean_doi=mean(duration),
              sd_doi=sd(duration),
              n=length(duration),
              CI_doi=1.96*sd(duration)/sqrt(n)) %>% 
    ggplot(aes(infection_id,mean_doi,color=case))+
    geom_line()+
    geom_errorbar(aes(x=infection_id,ymin=mean_doi-CI_doi,ymax=mean_doi+CI_doi))+
    # geom_smooth(method='nls',formula = y~exp(a + b * x), method.args=list(start = list(a = 0, b = 0)),se=F)+
    scale_color_manual(values=c('black',col,'gray'))+
    scale_y_continuous(limits = c(0,360), breaks = seq(0,360,60))+
    manuscript_theme+theme(axis.title = element_blank())
  return(panel)
}
  
panel_A <- make_panel_fig4('06','S',scenario_cols[1])
panel_B <- make_panel_fig4('06','G',scenario_cols[2])
panel_C <- make_panel_fig4('06','N',scenario_cols[3])
panel_D <- epi_results %>% 
  filter(PS=='06') %>% 
  filter(case=='R') %>% 
  mutate(scenario=factor(scenario, levels=c('S','G','N'))) %>% 
  group_by(scenario,infection_id) %>% 
  summarise(mean_doi=mean(duration),
            sd_doi=sd(duration),
            n=length(duration),
            CI_doi=1.96*sd(duration)/sqrt(n)) %>% 
  ggplot()+
  geom_line(aes(infection_id,mean_doi,color=scenario))+
  geom_errorbar(aes(x=infection_id,ymin=mean_doi-CI_doi,ymax=mean_doi+CI_doi,color=scenario))+
  # geom_smooth(aes(infection_id,mean_doi,color=scenario), method='nls',formula = y~exp(a + b * x), method.args=list(start = list(a = 0, b = 0)),se=F)+
  # geom_smooth(aes(infection_id,mean_doi-CI_doi,color=scenario), method='nls',formula = y~exp(a + b * x), method.args=list(start = list(a = 0, b = 0)),se=F, alpha=0.5, linetype='dashed',size=0.5)+
  # geom_smooth(aes(infection_id,mean_doi+CI_doi,color=scenario), method='nls',formula = y~exp(a + b * x), method.args=list(start = list(a = 0, b = 0)),se=F, alpha=0.5, linetype='dashed',size=0.5)+
  scale_color_manual(values=scenario_cols)+
  scale_x_continuous(limits=c(0,min(infection_limit$infection_limit)))+
  scale_y_continuous(limits = c(0,360), breaks = seq(0,360,60))+
  manuscript_theme+theme(axis.title = element_blank())
Fig <- plot_grid(panel_A,panel_B,panel_C,panel_D, labels=c('A','B','C','D'), ncol=2, align='vh', label_size = 18, scale=0.95)
y.grob <- textGrob("Mean duration of infection", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
x.grob <- textGrob("Number of infections", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
Fig_4 <- grid.arrange(arrangeGrob(Fig, left = y.grob, bottom = x.grob))


dev.off()
png('/home/shai/Dropbox/PLoS Biol/Fig_DOI_high.png', 4480*1.5,4490*1.5,units = 'px', res = 600)
# pdf('/home/shai/Dropbox/PLoS Biol/Fig_DOI_high.pdf', 16,12)
grid.arrange(arrangeGrob(Fig, left = y.grob, bottom = x.grob))
dev.off()


# Fit epidemiological curves ----------------------------------------------

# To compare between scenarios we fit the random case in which infection is
# regardless of modules to the first 50 infections. This tests the hypothesis
# that the decline in duration of infection is different between the scenarios.
# It does not test the effect of within vs between module infections on doi. Yet
# this last effect is irrelevant because in nature infections are not limited to
# either within or between modules. The first step is to test if an exponential
# model fits better than a linear one. It does. So the next step is to compare
# the exponential fits of the three scenarios. We do that by comparing the b
# coefficient, which is the one affectin gth decline.

library(minpack.lm)
epi_results

### Test several models to fit the data
model_comparison_results <- list()
for (scen in c('S','G','N')){
  # Model with all data ignoring group
  fit <- nlsLM(formula = duration~a*exp(-b*infection_id), 
               data = epi_results, 
               subset = (scenario==scen & case=='R' & infection_id<=50),
               start = c(a=1, b=0))
  # Linear model to show the decay is probably exponential
  fit_lm <- lm(duration~infection_id,
                   data=epi_results,
                   subset = (scenario==scen & case=='R' & infection_id<=50))
  # Compare models
  AIC <- AIC(fit,fit_lm)
  AIC <- AIC[order(AIC$AIC),]
  AIC$delta <- c(0, diff(AIC$AIC))
  AIC
  model_comparison_results[[scen]] <- AIC
}
model_comparison_results

# Explore the best model for each scenario
fit_best_model_scenario <- function(scen){
  fit <- nlsLM(formula = duration~a*exp(-b*infection_id), 
               data = epi_results, 
               subset = (scenario==scen & case=='R' & infection_id<=50),
               start = c(a=1, b=0))
  fit_summary <- summary(fit)
  a_coeff <- fit_summary$coefficients[1,1]
  b_coeff <- fit_summary$coefficients[2,1]
  b_coeff_se <- fit_summary$coefficients[2,2]
  return(c(a_coeff,b_coeff,b_coeff_se))
}

b_coeff_S <- fit_best_model_scenario('S')[2]
b_coeff_G <- fit_best_model_scenario('G')[2]
b_coeff_N <- fit_best_model_scenario('N')[2]

b_coeff_G/b_coeff_S
b_coeff_N/b_coeff_S


a_coeff_S <- fit_best_model_scenario('S')[1]
a_coeff_G <- fit_best_model_scenario('G')[1]
a_coeff_N <- fit_best_model_scenario('N')[1]

b_coeff_S_se <- fit_best_model_scenario('S')[3]
b_coeff_G_se <- fit_best_model_scenario('G')[3]
b_coeff_N_se <- fit_best_model_scenario('N')[3]


# Calculate the doi of a new infection in a 5-months old baby. The number of
# infections in 5 months depends on the EIR.
EIR_df <- data.frame(scenario='',run=0,EIR_run=0,stringsAsFactors = F)
for (i in 1:50){
  tmp <- get_data(parameter_space = '06', scenario = 'S', experiment = '001', run = i, cutoff_prob = 0.85, use_sqlite = T, tables_to_get = 'summary_general')[[1]]
  EIR_df <- rbind(EIR_df, c('S',i,mean(tmp$EIR)))
  tmp <- get_data(parameter_space = '06', scenario = 'G', experiment = '001', run = i, cutoff_prob = 0.85, use_sqlite = T, tables_to_get = 'summary_general')[[1]]
  EIR_df <- rbind(EIR_df, c('G',i,mean(tmp$EIR)))
  tmp <- get_data(parameter_space = '06', scenario = 'N', experiment = '001', run = i, cutoff_prob = 0.85, use_sqlite = T, tables_to_get = 'summary_general')[[1]]
  EIR_df <- rbind(EIR_df, c('N',i,mean(tmp$EIR)))
}
EIR_df <- EIR_df[-1,]
EIR_df$EIR_run <- as.numeric(EIR_df$EIR_run)
EIR_df <- as.tibble(EIR_df) %>% group_by(scenario) %>% summarise(EIR_mean=mean(EIR_run))

# Calculate per month
doi_by_month_S <- tibble(m=1:60, scenario='S')
doi_by_month_S$d <- a_coeff_S*exp(-b_coeff_S*subset(EIR_df, scenario=='S')$EIR_mean*doi_by_month_S$m)
doi_by_month_S$d_se_U <- a_coeff_S*exp(-(b_coeff_S+b_coeff_S_se)*subset(EIR_df, scenario=='S')$EIR_mean*doi_by_month_S$m)
doi_by_month_S$d_se_L <- a_coeff_S*exp(-(b_coeff_S-b_coeff_S_se)*subset(EIR_df, scenario=='S')$EIR_mean*doi_by_month_S$m)
doi_by_month_G <- tibble(m=1:60, scenario='G')
doi_by_month_G$d <- a_coeff_G*exp(-b_coeff_G*subset(EIR_df, scenario=='G')$EIR_mean*doi_by_month_G$m)
doi_by_month_G$d_se_U <- a_coeff_G*exp(-(b_coeff_G+b_coeff_G_se)*subset(EIR_df, scenario=='G')$EIR_mean*doi_by_month_S$m)
doi_by_month_G$d_se_L <- a_coeff_G*exp(-(b_coeff_G-b_coeff_G_se)*subset(EIR_df, scenario=='G')$EIR_mean*doi_by_month_S$m)
doi_by_month_N <- tibble(m=1:60, scenario='N')
doi_by_month_N$d <- a_coeff_N*exp(-b_coeff_N*subset(EIR_df, scenario=='N')$EIR_mean*doi_by_month_N$m)
doi_by_month_N$d_se_U <- a_coeff_N*exp(-(b_coeff_N+b_coeff_N_se)*subset(EIR_df, scenario=='N')$EIR_mean*doi_by_month_N$m)
doi_by_month_N$d_se_L <- a_coeff_N*exp(-(b_coeff_N-b_coeff_N_se)*subset(EIR_df, scenario=='N')$EIR_mean*doi_by_month_N$m)

doi_by_month <- rbind(doi_by_month_S,doi_by_month_G,doi_by_month_N)
doi_by_month %>% filter(m==5)
doi_by_month %>% 
  # filter(scenario=='G') %>% 
ggplot(aes(color=scenario))+
  geom_line(aes(x=m, y=d))+
  geom_line(aes(x=m, y=d_se_U), linetype='dashed')+
  geom_line(aes(x=m, y=d_se_L), linetype='dashed')+
  scale_color_manual(values = c('blue','orange','red'))+
  scale_x_continuous(breaks=seq(0,60,5))+
  scale_y_continuous(breaks=seq(0,300,50))+
  labs(x='Host age (months)', y='Predicted duration of infection')+
  manuscript_theme




make_panel_allele_curve <- function(ps,scen,col){
  panel <- epi_results %>% 
    filter(PS==ps) %>% 
    filter(scenario==scen) %>% 
    group_by(case,infection_id) %>% 
    summarise(mean_alleles=mean(alleles_seen),
              sd_alleles=sd(alleles_seen),
              n=length(alleles_seen),
              CI_alleles=1.96*sd(alleles_seen)/sqrt(n)) %>% 
    ggplot(aes(infection_id,mean_alleles,color=case))+
    geom_line()+
    geom_errorbar(aes(x=infection_id,ymin=mean_alleles-CI_alleles,ymax=mean_alleles+CI_alleles))+
    # geom_smooth(method='nls',formula = y~exp(a + b * x), method.args=list(start = list(a = 0, b = 0)),se=F)+
    scale_color_manual(values=c('black',col,'gray'))+
    # scale_y_continuous(limits = c(0,360), breaks = seq(0,360,60))+
    manuscript_theme+theme(axis.title = element_blank())
  return(panel)
}

panel_A <- make_panel_allele_curve('06','S',scenario_cols[1])
panel_B <- make_panel_allele_curve('06','G',scenario_cols[2])
panel_C <- make_panel_allele_curve('06','N',scenario_cols[3])
panel_D <- epi_results %>% 
  filter(PS=='06') %>% 
  filter(case=='R') %>% 
  mutate(scenario=factor(scenario, levels=c('S','G','N'))) %>% 
  group_by(scenario,infection_id) %>% 
  summarise(mean_alleles=mean(alleles_seen),
            sd_alleles=sd(alleles_seen),
            n=length(alleles_seen),
            CI_alleles=1.96*sd(alleles_seen)/sqrt(n)) %>% 
  ggplot()+
  geom_line(aes(infection_id,mean_alleles,color=scenario))+
  geom_errorbar(aes(x=infection_id,ymin=mean_alleles-CI_alleles,ymax=mean_alleles+CI_alleles,color=scenario))+
  # geom_smooth(aes(infection_id,mean_doi,color=scenario), method='nls',formula = y~exp(a + b * x), method.args=list(start = list(a = 0, b = 0)),se=F)+
  # geom_smooth(aes(infection_id,mean_doi-CI_doi,color=scenario), method='nls',formula = y~exp(a + b * x), method.args=list(start = list(a = 0, b = 0)),se=F, alpha=0.5, linetype='dashed',size=0.5)+
  # geom_smooth(aes(infection_id,mean_doi+CI_doi,color=scenario), method='nls',formula = y~exp(a + b * x), method.args=list(start = list(a = 0, b = 0)),se=F, alpha=0.5, linetype='dashed',size=0.5)+
  scale_color_manual(values=scenario_cols)+
  # scale_x_continuous(limits=c(0,min(infection_limit$infection_limit)))+
  # scale_y_continuous(limits = c(0,360), breaks = seq(0,360,60))+
  manuscript_theme+theme(axis.title = element_blank())
Fig <- plot_grid(panel_A,panel_B,panel_C,panel_D, labels=c('A','B','C','D'), ncol=2, align='vh', label_size = 18, scale=0.95)
y.grob <- textGrob("Mean # alleles accumulated", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
x.grob <- textGrob("Number of infections", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
grid.arrange(arrangeGrob(Fig, left = y.grob, bottom = x.grob))


# fit_scearios <- nlsLM(formula = duration~
#                         as.numeric(scenario=='S')*a1*exp(-b1*infection_id)+ # The as.numeroc part is a dummy variable to make this categorical
#                         as.numeric(scenario=='G')*a2*exp(-b2*infection_id)+ # The as.numeroc part is a dummy variable to make this categorical
#                         as.numeric(scenario=='N')*a3*exp(-b3*infection_id), # The as.numeroc part is a dummy variable to make this categorical
#                       data = epi_results,
#                       subset = (infection_id<=min(infection_limit$infection_limit) & case=='R'),
#                       start = c(a1=1,b1=0, a2=1,b2=0, a3=1,b3=0))
# # Linear model to show the decay is probably exponential
# fit_lm_scearios <- lm(duration~infection_id+scenario,
#                       data=epi_results,
#                       subset = (infection_id<=min(infection_limit$infection_limit)) & case=='R')
# 
# AIC(fit_scearios,fit_lm_scearios)
# summary(fit_scearios)
