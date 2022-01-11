#------------
# plot regression result
# Fig. 2 
# By Yang B. and Lin Y. 
# updated October 2021
#------------
######################################################
## daily_ct_rt: daily case counts/sample counts, incidence-based Rt; 
#               Ct-based Rt estimated from main model (with prediction interval)
##              from "4_model_loglinear"
######################################################
#
# load packages
require(ggplot2)
library(tidyverse)
library(pROC)
library(gridExtra)
require(e1071)
#
#setwd("/Users/vanialam/OneDrive - connect.hku.hk/vanialam/research_vania/epi_wave_2021/program/2021_09_R1/publish (EDIT HERE)/")
# read in "daily_ct_rt.csv"
ct.rt <- read.csv("results/daily_ct_rt.csv",as.is = T)
## variable explanations (for Ct-based Rt in the data)
# fit - point estimates 
# upr - upper range for prediction interval 
# lwr - lower range for prediction interval
#

source("ggplot_default.R")
# functions set for plotting 
predPlot = function(ct, period, panel){
        
        y.max = 6
        
        if(period == 'training') {
                data_used = ct %>%
                        filter(training == 1)
                dates = range(data_used$date)
        }
        
        if(period == 'testing') {
                data_used = ct %>%
                        filter(testing == 1)
                dates = range(data_used$date)
        }
        
        
        data_rt = data_used %>% # to truncate incidence-based Rt
                filter(date <= max(dates)-7)
        
        p = ggplot(data = data_used %>%
                           mutate(upr = ifelse(upr > y.max, y.max, upr))) +
                
                #epi curve
                geom_bar(data=data_used,
                         aes(x = date,
                             y = all.cases),
                         stat = 'identity',
                         fill = '#dad5d4') +
                
                #case rt
                geom_polygon(data = tibble(
                        date_new = c(data_rt$date, 
                                     rev(data_rt$date)),
                        y_new = c(data_rt$local.rt.lower,
                                  rev(data_rt$local.rt.upper))) %>%
                                mutate(y_new = ifelse(y_new > y.max, y.max, y_new)),
                        aes(x = date_new,
                            y = y_new*25,
                            fill = 'empirical'),
                        color = NA,
                        alpha = 0.2) +
                
                geom_line(data = data_rt ,
                          aes(x = date,
                              y = local.rt.mean*25,
                              color = 'empirical'),
                          size = 1) +
                #ct rt
                geom_point(aes(x = date,
                               y = fit*25,
                               color = 'predicted'),
                           size = 2) +
                
                geom_segment(aes(x = date,
                                 y = lwr*25,
                                 xend = date,
                                 yend = upr*25,
                                 color = 'predicted'),
                             size = 0.8) + 
                
                scale_y_continuous(name = 'Cases',
                                   limits = c(0, 150),
                                   expand = c(0, 0),
                                   breaks = seq(0, 150, 30),
                                   position = 'right',
                                   sec.axis = sec_axis(~./25, 
                                                       name = 'Rt')) +
                
                scale_x_date(name = 'Date',
                             limits = c(dates[1]-1,dates[2]+1),
                             date_breaks = "1 week", 
                             date_labels = "%d/%m",
                             expand = c(0.01, 0.01))  +
                
                geom_hline(yintercept = 1*25,
                           linetype = 'dashed',
                           size = 1,
                           color = 'grey') +
                
                theme(legend.position = c(0.8, 0.85),
                      legend.title = element_text(size = 16, face = 'bold')) +
                labs(title = panel)
        
        if(period == 'training'){
                
                p = p +
                        scale_fill_manual(name = NULL,
                                          values = c(empirical = 'black',
                                                     predicted = '#9292e4'),
                                          labels = c('Incidence-based', 'Ct-based, training'), 
                                          guide = "none") +
                        scale_color_manual(name = NULL,
                                           values = c(empirical = 'black',
                                                      predicted = '#9292e4'),
                                           labels = c('Incidence-based', 'Ct-based, training'))
                
        }
        
        if(period == 'testing'){
                
                p = p +
                        scale_fill_manual(name = NULL,
                                          values = c(empirical = 'black',
                                                     predicted = '#e49292'),
                                          labels = c('Incidence-based', 'Ct-based, testing'), 
                                          guide = "none") +
                        scale_color_manual(name = NULL,
                                           values = c(empirical = 'black',
                                                      predicted = '#e49292'),
                                           labels = c('Incidence-based', 'Ct-based, testing'))
                
        }
        
        p
        
}
#
predBoxPlot = function(ct, panel){
        
        data_used = ct %>%
                filter(training == 1 | testing == 1) %>%
                filter(!is.na(fit)) %>%
                mutate(
                        pred_rt_cat = factor(cut(fit,
                                                 c(-1, 0.5, 1, 1.5, 10))),
                        group = factor(ifelse(training == 1,
                                              'Training, wave 3',
                                              'Testing, wave 4'),
                                       levels = c('Training, wave 3',
                                                  'Testing, wave 4'))
                )
        
        ggplot(data = data_used) +
                geom_boxplot(aes(
                        x = pred_rt_cat,
                        y = local.rt.mean,
                        fill = group),
                        width = 0.3) +
                geom_hline(yintercept = 1,
                           linetype = 'dashed',
                           size = 1,
                           color = 'grey') +
                scale_y_continuous(name = 'Incidence-based Rt',
                                   limits = c(0, 4),
                                   breaks = seq(0, 4, 1),
                                   expand = c(0, 0)) +
                scale_x_discrete(name = 'Ct-based Rt',
                                 expand = c(0.01, 0.01),
                                 labels = c('<0.5', '0.5-1.0', '1.0-1.5', '>1.5')) +
                scale_fill_manual(name = '',
                                  values = c('#c9c9ff', '#f3cfcf')) +
                theme(legend.position = c(0.35, 0.9)) +
                labs(title = panel)
}
#
retroExamplePlot = function(df, date, panel){
        
        date_max = date + 6
        date_cut = date - 1
        date_0 = as.Date('2020-11-19')
        date_end = as.Date('2021-02-01')
        
        y.max = 6
        
        data_used = df %>%
                filter(date >=  date_0 &
                               date <= date_cut)
        
        data_ct = df %>%
                filter(date >=  date_0 &
                               date <= date_max)
        
        data_case = df %>%
                filter(date >=  date_0 &
                               date <= date_max)
        
        
        p = ggplot() +
                
                # epi curve
                geom_bar(data=data_case,
                         aes(x = date,
                             y = all.cases),
                         stat = 'identity',
                         fill = '#dad5d4') +
                
                geom_hline(yintercept = 1 * 20,
                           linetype = 'dashed',
                           size = .8,
                           color = 'black') +
                
                # case rt
                geom_polygon(data = tibble(
                        date_new = c(data_used$date, 
                                     rev(data_used$date)),
                        y_new = c(data_used$local.rt.lower,
                                  rev(data_used$local.rt.upper))) %>%
                                mutate(y_new = ifelse(y_new > y.max, y.max, y_new)),
                        aes(x = date_new,
                            y = y_new * 20,
                            fill = 'empirical'),
                        color = NA,
                        alpha = 0.2) +
                
                geom_line(data = data_used,
                          aes(x = date,
                              y = local.rt.mean * 20,
                              color = 'empirical'),
                          size = .8) +
                # ct rt
                geom_point(data = data_ct,
                           aes(x = date,
                               y = fit * 20,
                               color = 'predicted'),
                           size = 1.3,
                           alpha = 0.8) +
                
                geom_segment(data = data_ct,
                             aes(x = date,
                                 y = lwr * 20,
                                 xend = date,
                                 yend = upr * 20,
                                 color = 'predicted'),
                             size = .7,
                             alpha = 0.8) +
                
                ### arrow indicating current day for prediction
                geom_segment(aes(x = max(data_ct$date), y = 113, 
                                 xend = max(data_ct$date), yend = 100),
                             arrow = arrow(length = unit(0.25, "cm")),
                             size = .8) +
                
                scale_y_continuous(name = 'Cases',
                                   limits = c(0, 120),
                                   expand = c(0, 0),
                                   breaks = seq(0, 120, 20),
                                   position = 'right',
                                   sec.axis = sec_axis(~./20, 
                                                       name = 'Rt')) +
                
                scale_x_date(name = 'Date of report',
                             limits = c(date_0-0.5, date_end+0.5),
                             date_breaks = "1 week", 
                             date_labels = "%d/%m",
                             expand = c(0.01, 0)) +
                
                theme(legend.title = element_text(size = 16, face = 'bold'),
                      plot.title.position = 'plot',
                      plot.title = element_text(size = 16, face = 'bold')) +
                labs(title = panel) +
                annotate("text", x = max(data_ct$date)-6.5, y = 117, 
                         label = "Date of estimation", size=4) +
                scale_fill_manual(name = NULL,
                                  values = c(empirical = 'black',
                                             predicted = '#e49292'),
                                  labels = c('Incidence-based', 'Ct-based, testing'), 
                                  guide = "none") +
                scale_color_manual(name = NULL,
                                   values = c(empirical = 'black',
                                              predicted = '#e49292'),
                                   labels = c('Incidence-based', 'Ct-based, testing'))
        
        if(panel == 'a'){
                
                p = p + theme(legend.position = c(0.65, 0.85))
                
        } else {
                
                p = p + theme(legend.position = 'none')
                
        }
        p
}
###
#-------
date_max <- as.Date("2021-03-31")
ct.rt$date <- as.Date(ct.rt$date)
ct.rt2 <- ct.rt%>%
        
        mutate(training = ifelse(date >= as.Date('2020-07-01') &
                                         date <= as.Date('2020-08-31'),
                                 1,0),
               testing = ifelse(date >= as.Date('2020-11-01') &
                                        date <= as.Date('2021-03-31'),
                                          1,0) 
               
              ) %>% 
        filter(date >= as.Date('2020-07-01') &
                             
                       date <= date_max
               )
#
ct.rt2$upr[ct.rt2$upr>6&!is.na(ct.rt2$upr)] <- 6 # y-limit for panel a is 6

p = grid.arrange(
        grobs = list(retroExamplePlot(ct.rt2, date = as.Date('2020-12-01'), 'a'),
                     retroExamplePlot(ct.rt2, date = as.Date('2020-12-15'), ''),
                     retroExamplePlot(ct.rt2, date = as.Date('2021-01-10'), ''),
                     retroExamplePlot(ct.rt2, date = as.Date('2021-01-25'), ''),
                     predPlot(ct.rt2, 'training', 'b'),
                     predPlot(ct.rt2, 'testing', 'c'),
                     predBoxPlot(ct.rt2, 'd')),
        heights = c(1, 1, 1, 1),
        widths = c(5, 7, 4),
        layout_matrix = rbind(c(1,5,7),
                              c(2,5,7),
                              c(3,6,6),
                              c(4,6,6))
)
ggsave("results/Fig_2.pdf",p,width = 22, height = 12)
##
######

## end of script

#####