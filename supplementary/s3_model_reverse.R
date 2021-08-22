#------------
# reverse validation
# ED.Fig 4 in paper
#------------
#
# load packages
require(ggplot2)
library(tidyverse)
library(pROC)
library(gridExtra)
require(e1071)
#
######################################################
## data_daily_all: daily case counts/sample counts, incidence-based Rt; 
##                 daily Ct mean, median and skewness (imputed)
#######################################################
#
# read in "data_daily_all.csv"
daily.linelist <- read.csv("data_daily_all.csv",as.is=T)
#
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
        
        p = ggplot(data = data_used %>%
                           mutate(upr = ifelse(upr > y.max, y.max, upr))) +
                
                geom_polygon(data = tibble(
                        date_new = c(data_used$date, 
                                     rev(data_used$date)),
                        y_new = c(data_used$local.rt.lower,
                                  rev(data_used$local.rt.upper))) %>%
                                mutate(y_new = ifelse(y_new > y.max, y.max, y_new)),
                        aes(x = date_new,
                            y = y_new,
                            fill = 'empirical'),
                        color = NA,
                        alpha = 0.2) +
                
                geom_line(aes(x = date,
                              y = local.rt.mean,
                              color = 'empirical'),
                          size = 1) +
                
                geom_point(aes(x = date,
                               y = fit,
                               color = 'predicted'),
                           size = 2) +
                
                geom_segment(aes(x = date,
                                 y = lwr,
                                 xend = date,
                                 yend = upr,
                                 color = 'predicted'),
                             size = 0.8) +
                
                scale_y_continuous(name = 'Rt',
                                   limits = c(0, y.max),
                                   expand = c(0, 0),
                                   breaks = seq(0, y.max, 1)) +
                scale_x_date(name = 'Date',
                             limits = dates,
                             date_breaks = "1 week", 
                             date_labels = "%d/%m",
                             expand = c(0.01, 0.01))  +
                
                geom_hline(yintercept = 1,
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
                                              'Training, wave 4',
                                              'Testing, wave 3'),
                                       levels = c('Training, wave 4',
                                                  'Testing, wave 3'))
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
#-------
date_max <- as.Date("2021-03-31")
daily.linelist$date <- as.Date(daily.linelist$date)

ct.rt2 = daily.linelist %>%
        mutate(training = ifelse(date >= as.Date('2020-11-20') &
                                         date <= as.Date('2020-12-31'),
                                 1,
                                 0),
               testing = ifelse(date >= as.Date('2020-07-01') &
                                        date <= as.Date('2020-08-31'),
                                1,
                                0)
        ) %>%
        filter(date >= as.Date('2020-07-01') &
                       date <= date_max)

m = lm(log(local.rt.mean) ~ mean + skewness.imputed, 
       data = ct.rt2[ct.rt2$training == 1, ])

pred = predict(m, newdata = ct.rt2, interval = 'prediction')

ct.rt2 = cbind(ct.rt2, exp(pred))


p = grid.arrange(
        grobs = list(
                predPlot(ct.rt2, 'training', 'a'),
                predBoxPlot(ct.rt2, 'c'),
                predPlot(ct.rt2, 'testing', 'b')),
        widths = c(2, 1),
        heights = c(1, 1),
        layout_matrix = rbind(c(1, 2),
                              c(3, 3))
)
## export results 
ggsave('ED_Fig_4.pdf',p,width = 12, height = 8)
##
#####

## end of script

#####
