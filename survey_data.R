library(tidyverse)
library(cowplot)

mytheme_no_legend <- theme_bw() + theme(
  legend.title  = element_text(colour = "black", size=17),
  legend.position = "none",
  legend.key = element_blank(),
  legend.text  = element_text(colour = "black", size=17),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  panel.border = element_blank(),
  # panel.border = element_rect(colour = "black", size=1.3),
  axis.text = element_text(color='black', family="Helvetica", size=14),
  axis.title = element_text(color='black', family="Helvetica", size=14),
  strip.text.x = element_text(family = "Helvetica", size = 14),
  strip.text.y = element_text(family = "Helvetica", size = 14),
  axis.ticks = element_line(size = 1.3),
  strip.background = element_rect( fill = "transparent", size = 1.3, colour = "black"  ),
  strip.text = element_text(size = 19)
)
survey_data <- tibble(ID=1:8,
                      survey=c('S1',
                               'S2',
                               'IRS1',
                               'S3/IRS2',
                               'S4',
                               'IRS3',
                               'S5',
                               'S6'),
                      Date=c('Oct. 2012',
                             'May-Jun. 2013',
                             'Oct-Dec 2013',
                             'May-Jun. 2014',
                             'Oct. 2014',
                             'Jan-Febuary 2015',
                             'Oct. 2015',
                             'Jun. 2016'),
                      season=c('EOW','EOD','EOW','EOD','EOW','MOD','EOW','EOD'),
                      Prevalence=c(808/1923,513/1902,NA,535/1822,430/1866,NA,545/2022,272/2091),
                      var_types=c(35345,27586,NA,24962,16222,NA,19630,18103))

Fig_prev <- ggplot(survey_data, aes(x=ID,y=Prevalence))+
  geom_rect(aes(xmin = 0.5, xmax = 1.5,ymin = -Inf, ymax = Inf),fill = "dark green", alpha=0.1)+
  geom_rect(aes(xmin = 1.5, xmax = 2.5,ymin = -Inf, ymax = Inf),fill = "brown", alpha=0.1)+
  geom_rect(aes(xmin = 2.5, xmax = 3.5,ymin = -Inf, ymax = Inf),fill = "dark green", alpha=0.1)+
  geom_rect(aes(xmin = 3.5, xmax = 4.5,ymin = -Inf, ymax = Inf),fill = "brown", alpha=0.1)+
  geom_rect(aes(xmin = 4.5, xmax = 5.5,ymin = -Inf, ymax = Inf),fill = "dark green", alpha=0.1)+
  geom_rect(aes(xmin = 5.5, xmax = 6.5,ymin = -Inf, ymax = Inf),fill = "brown", alpha=0.1)+
  geom_rect(aes(xmin = 6.5, xmax = 7.5,ymin = -Inf, ymax = Inf),fill = "dark green", alpha=0.1)+
  geom_rect(aes(xmin = 7.5, xmax = 8.5,ymin = -Inf, ymax = Inf),fill = "brown", alpha=0.1)+
  geom_point(size=5)+
  geom_line()+
  geom_vline(xintercept = c(3,4,6), linetype='dashed')+
  scale_x_continuous(breaks=survey_data$ID, labels = survey_data$survey)+
  scale_y_continuous(breaks=seq(0,0.5,0.1),labels=str_pad(seq(0.,0.5,0.1),4,'right',0))+
  mytheme_no_legend+theme(axis.title.x = element_blank())
  
Fig_var <- ggplot(survey_data, aes(x=ID,y=var_types))+
  geom_rect(aes(xmin = 0.5, xmax = 1.5,ymin = -Inf, ymax = Inf),fill = "dark green", alpha=0.1)+
  geom_rect(aes(xmin = 1.5, xmax = 2.5,ymin = -Inf, ymax = Inf),fill = "brown", alpha=0.1)+
  geom_rect(aes(xmin = 2.5, xmax = 3.5,ymin = -Inf, ymax = Inf),fill = "dark green", alpha=0.1)+
  geom_rect(aes(xmin = 3.5, xmax = 4.5,ymin = -Inf, ymax = Inf),fill = "brown", alpha=0.1)+
  geom_rect(aes(xmin = 4.5, xmax = 5.5,ymin = -Inf, ymax = Inf),fill = "dark green", alpha=0.1)+
  geom_rect(aes(xmin = 5.5, xmax = 6.5,ymin = -Inf, ymax = Inf),fill = "brown", alpha=0.1)+
  geom_rect(aes(xmin = 6.5, xmax = 7.5,ymin = -Inf, ymax = Inf),fill = "dark green", alpha=0.1)+
  geom_rect(aes(xmin = 7.5, xmax = 8.5,ymin = -Inf, ymax = Inf),fill = "brown", alpha=0.1)+
  geom_point(size=5, color='purple')+
  geom_line(color='purple')+
  geom_vline(xintercept = c(3,4,6), linetype='dashed')+
  scale_x_continuous(breaks=survey_data$ID, labels = survey_data$survey)+
  scale_y_continuous(breaks=seq(15000,40000,5000),labels=seq(15000,40000,5000))+
  labs(y='var DBLa types sampled')+
  mytheme_no_legend+theme(axis.title.x = element_blank())

Fig <- plot_grid(Fig_prev, Fig_var, nrow=2, align='vh')

png('~/Dropbox/Images for presentations/BGU Job talk/survey_data.png', 1600,1200,res=250)
Fig
dev.off()



