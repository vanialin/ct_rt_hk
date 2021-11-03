#-----------
# sourced functions
#-----------

##
### source function for justifying the best training period ----
# this function return the most suitable time period 
# for training the regression model
select_training_period <- function(df){
        initial_start <- df$date[which(df$count>15)][1]
        interval <- 2:6*10 
        r.sqr.mat <- matrix(NA,60,5)
        for (i in 1:60){ # 60 starts
               for (j in 1:length(interval)){
                       start_date <- initial_start+i
                       end_date <- start_date+interval[j]-1
                       training.tmp <- 
                               df[df$date%in%seq(start_date,end_date,1),]
                       ##
                       lm.tmp <- lm(log(rt_est)~mean+skewness_imputed,
                                    data=training.tmp)
                       
                       r.sqr.mat[i,j] <- summary(lm.tmp)$adj.r.square
               }
        }
        
        ## select based on adjusted R-square
        fit <- which(r.sqr.mat == max(r.sqr.mat), arr.ind = TRUE)
        final_training_period <- 
                seq(initial_start+fit[1],
                    initial_start+fit[1]+interval[fit[2]],1)
        df.train <- df[df$date%in%final_training_period,]
        
        ###
        vec.tmp <- c(NA,9)
        cor.test.mean <- with(df.train,cor.test(mean,log(rt_est),
                                                method = "spearman",
                                                use="na.or.complete"))
        vec.tmp[1:2] <- round(c(cor.test.mean$estimate,cor.test.mean$p.value),3)
        cor.test.skewness <- with(df.train,cor.test(skewness,log(rt_est),
                                                    method = "spearman",
                                                    use="na.or.complete"))
        vec.tmp[3:4] <- round(c(cor.test.skewness$estimate,
                                cor.test.skewness$p.value),3)
        lm.final <- lm(log(rt_est)~mean+skewness_imputed,data=df.train)
        vec.tmp[5:8] <- c(summary(lm.final)$coefficients[2,c(1,4)],
                          summary(lm.final)$coefficients[3,c(1,4)])
        vec.tmp[9] <- summary(lm.final)$adj.r.square 
        vec.tmp[c(5,7)] <- exp(vec.tmp[c(5,7)])
        names(vec.tmp) <- c("rho_mean","p_mean","rho_skewness","p_skewness",
                            "coef_mean","p_coef_mean",
                            "coef_skewness","p_coef_skewness",
                            "adj_r_sqr")
        
        vec.tmp <- round(vec.tmp,3)
        
        return(list(period=final_training_period[c(1,length(final_training_period))],
                    characteristics=vec.tmp))
}

#
#

### source function for evaluating consistency between Rt ----
Rt_consistency <- function(df,rt1,rt2){
        # rt1 - incidence-based Rt
        # rt2 - Ct-based Rt
        var_name <- c(rt1,rt2)
        cor.rt <- matrix(NA,2,5)
        for (i in 1:2){
                cor.rt[i,1:3] <- 
                        SpearmanRho(log(df[,var_name[i]]),log(df$Rt),
                                    use="na.or.complete",conf.level = .95)
                cor.rt[i,4] <- sum(abs(log(df[,var_name[i]])-log(df$Rt)))/nrow(df)
        }
        cor.rt[2,5] <- 
                round((nrow(df[((df[,rt2]-1)*(df$Rt-1))>0,])/
                               nrow(df))*100,1)
        
        colnames(cor.rt) <- c("rho","lower.ci","uppper.ci","mae","prop.")
        return(cor.rt)
}

#
#

### for plotting epi curve with Rt -----
epi_plot <- function(df,truth=F,training=F){
        y_upper <- round(max(df$count)+50,digits = -2)
        times <- y_upper/3
        
        p <- ggplot(df,aes(x=date)) + 
                # epi curve
                geom_bar(aes(y = count),
                         stat = 'identity',
                         fill = '#dad5d4') +
                # reference Rt at 1
                geom_hline(yintercept = 1 * times,
                           linetype = 'dashed',
                           size = .8,
                           color = 'black') +
                # case rt
                geom_ribbon(aes(ymin = rt_lb * times,
                                ymax = rt_ub * times, fill = "red"),
                            color = NA,
                            alpha = 0.2,
                            show.legend = F) +
                geom_line(aes(y = rt_est * times),
                          color = "red",
                          alpha = 0.6) +
                scale_y_continuous(name = 'Cases',
                                   limits = c(0, y_upper),
                                   expand = c(0, 0),
                                   breaks = seq(0, y_upper, length.out = 6),
                                   position = 'right',
                                   sec.axis = sec_axis(~./times, 
                                                       name = 'Rt'))
        
        if (truth == T){
                p <- p + geom_line(aes(y = Rt * times),
                                   color = "black", alpha = 0.6)
        }
        
        ##
        if (training==T){
                p <- p + 
                        geom_vline(xintercept = max(df$date[df$period==1]),
                                   linetype = 'dashed',
                                   size = .6,
                                   color = 'grey') +
                        geom_bracket(xmin = min(df$date), 
                                     xmax = max(df$date[df$period==1]), 
                                     y.position = 2.5 * times,
                                     label = "Training",
                                     label.size = 4,
                                     vjust = -.3) +
                        theme(legend.title = element_text(size = 10))
        }       
        
        return(list(epi=p,p_limit=y_upper))        
        
}

#
#

### for merging Ct-Rt -----
merge_Ct_Rt <- function(rt0,ct0,smoothing){
        ct1 <- ct0 %>% filter(ct_obs < 40) %>% 
                group_by(sampled_time) %>% 
                summarise(count=n(),
                          mean=mean(ct_obs),
                          skewness=skewness(ct_obs)) %>%
                ungroup() 
        
        complete_days <- seq(0,max(ct1$sampled_time),1)
        dummy_row <- tibble(sampled_time=complete_days[!complete_days%in%ct1$sampled_time],
                            count=0,mean=NA,skewness=NA)
        ct2 <- bind_rows(ct1,dummy_row) %>% arrange(sampled_time)
        
        start_date <- as.Date("2020-07-01")
        ct3 <- ct2 %>% 
                mutate( ## impute skewness
                        skewness_imputed = sapply(ct2$sampled_time, function(d){
                                mean(ct2$skewness[ct2$sampled_time %in% c((d-7):(d-1))], 
                                     na.rm = TRUE)
                        }),
                        ## If calculating moving average by taking mean of ALL raw Cts
                        mean_7d=sapply(sampled_time,function(d){
                                mean(ct0$ct_obs[ct0$ct_obs<40&
                                                        ct0$sampled_time%in%((d-6):d)])
                        }),
                        mean_14d=sapply(sampled_time,function(d){
                                mean(ct0$ct_obs[ct0$ct_obs<40&
                                                        ct0$sampled_time%in%((d-13):d)])
                        }),
                        ## If calculating by taking mean of mean Ct over the time window
                        #mean_d7=zoo::rollmean(mean, 7,align="right",fill=NA),
                        #mean_d14=zoo::rollmean(mean, 14,align="right",fill=NA),
                        date=start_date+sampled_time
                ) 
        ct3$skewness_imputed[!is.na(ct3$skewness)] <- ct3$skewness[!is.na(ct3$skewness)]
        
        rt1 <- rt0 %>% mutate(date=as.Date(date)) %>%
                dplyr::select(date,mean,lower_95,upper_95) %>%
                rename(rt_est=mean,rt_lb=lower_95,rt_ub=upper_95)
        ct4 <- merge(ct3,rt1,by.x = "date",by.y = "date")
        ### out
        if (smoothing==T){
                return(ct4)
        } else {
                return(ct4=dplyr::select(ct4, -c("mean_7d", "mean_14d")))
        }
}

#
#

### for plotting incidence-based and Ct-based Rt ----
rt_plot <- function(df,time,add_epi){ 
        ###
        if (time == "all"){
                df_plot <- df 
        } else if (time == "testing"){
                df_plot <- df %>% filter(period == 2)
        }
        ##
        cols <- c(empirical = '#d9534f',
                  predicted = '#428bca',
                  truth = "#400000") 
        
        ### trying to add epidemic curves
        if (add_epi == T){
                y_upper <- round(max(df_plot$count)+50,digits = -2)
                times <- y_upper/3
                
                p_tmp <- df_plot %>% ggplot(aes(x=date)) +
                        geom_bar(aes(y = count),
                                 stat = 'identity',
                                 fill = '#dad5d4') +
                        ### reference line for Rt
                        geom_hline(yintercept = 1*times,
                                   linetype = 'dashed',
                                   size = .4,
                                   color = 'black')
        } else {
                p_tmp <- df_plot %>% ggplot(aes(x=date)) +
                        geom_hline(yintercept = 1,
                                   linetype = 'dashed',
                                   size = .4,
                                   color = 'black')
        }
        
        p_tmp <- p_tmp +
                
                ## simulation truth
                geom_line(aes(y = Rt*(times^(add_epi)),
                              color = 'truth'),
                          alpha = 0.8) +
                
                ## incidence-based Rt
                geom_ribbon(aes(ymin = rt_lb*(times^(add_epi)),
                                ymax = rt_ub*(times^(add_epi)), fill = 'empirical'),
                            color = NA,
                            alpha = 0.2,
                            show.legend = F) +
                geom_line(aes(y = rt_est*(times^(add_epi)),
                              color = 'empirical'),
                          alpha = 0.8) +
                
                ## Ct-Rt 
                geom_segment(aes(y = lwr*(times^(add_epi)),
                                 xend = date,
                                 yend = upr*(times^(add_epi)),
                                 color = 'predicted'),
                             size = .4,
                             alpha = 0.8) +
                geom_point(aes(y = fit*(times^(add_epi)),
                               color = 'predicted'),
                           alpha = 0.7,
                           size = .2) +
                
                ## legends
                scale_x_date(name = 'Date',
                             limits = c(min(df_plot$date)-0.5, max(df_plot$date)+0.5),
                             date_breaks = "1 month", 
                             date_labels = "%d/%m/%y",
                             expand = c(0.01, 0)) +
                scale_fill_manual(name = '',
                                  values = cols,
                                  labels = c('Incidence-Rt', 'Ct-Rt','Simulation truth'),
                                  guide = "none") +
                scale_color_manual(name = '',
                                   values = cols,
                                   labels = c('Incidence-Rt', 'Ct-Rt','Simulation truth')) + 
                theme(legend.position = c(0.85, 0.95))
        
        if (add_epi == T){
                p_tmp <- p_tmp + scale_y_continuous(name = 'Cases',
                                                    limits = c(0, y_upper),
                                                    expand = c(0, 0),
                                                    breaks = seq(0, y_upper, length.out = 6),
                                                    position = 'right',
                                                    sec.axis = sec_axis(~./times, 
                                                                        name = 'Rt'))
        } else {
                p_tmp <- p_tmp + 
                        scale_y_continuous(name = 'Rt',
                                           limits = c(0, 3.5),
                                           breaks = seq(0, 3, 1),
                                           expand = c(0, 0)) 
        }
        
        
        if (time == "all"){
                if (sum(df_plot$period==1)<40){
                        lab_for_training <- "Training"
                } else {
                        lab_for_training <- "Training period"
                }
                p_tmp <- p_tmp + 
                        geom_vline(xintercept = max(df_plot$date[df_plot$period==1]),
                                   linetype = 'dashed',
                                   size = .6,
                                   color = 'grey') +
                        geom_bracket(xmin = min(df_plot$date), 
                                     xmax = max(df_plot$date[df_plot$period==1]), 
                                     y.position = 2.8*(times^add_epi),
                                     label = lab_for_training,
                                     tip.length = 0.001)
        } else {
                p_tmp
        }
        
        return(p_tmp)
}

#
#

### function for plotting bars ----
for_ci_plot <- function(x,a,col){ 
        # x is horizontal location, a is vector defining point est (ci)
        lines(c(x-0.05,x+0.05),rep(a[2],2),col=col,lwd=1.3)
        lines(c(x-0.05,x+0.05),rep(a[3],2),col=col,lwd=1.3)
        lines(c(x-0.05,x+0.05),rep(a[1],2),col=col,lwd=2)
        lines(rep(x,2),a[2:3],col=alpha(col,.7),lwd=1.2)
}

##
#####

## end of script

#####