#-----------
# function to output 
# all model results
# with only version for daily
#-----------

evaluate_daily_funcs <- function(rt0,ct0,sim_dyn, # input data
                                 ## demanded by nested function
                                 extra, 
                                 add_legend){
        # sim_dyn = seir_dynamics (sourced, simulation truth)
        start_date <- as.Date("2020-07-01") # day 0 for simulation
        ct0$date <- start_date+ct0$sampled_time # to facilitate plotting
        #rt0
        # get summary df
        ct4 <- merge_Ct_Rt(rt0,ct0,smoothing = T)
        ## based on df (plot raw distribution)
        p1 <- ct4 %>%
                filter(date<=max(ct0$date)) %>% # to get consistent date scale
                epi_plot() 
        p1_out <- p1[[1]]
        
        p2 <- ct0 %>% filter(ct_obs<40) %>% ggplot() +
                geom_smooth(aes(date,ct_obs))+scale_y_reverse(name="Ct value")
        
        grid.newpage()
        ### OUTPUT 1: CORRELATION PLOT ###
        ## Create a graphical object, but not draw here
        p <- rbind(ggplotGrob(p1_out), ggplotGrob(p2), size = "first") 
        
        ##
        #### correlation and regression 
        # define training period
        result <- select_training_period(ct4,time_frame = "1d",
                                         extra = extra,include_skew = T)
        training_period <- seq(result$period[1],
                               result$period[2],1)
        
        ct4$period <- ifelse(ct4$date%in%training_period,1,
                             # 1-training; 2-testing
                             ifelse(ct4$date<min(training_period),0,2))
        
        correlation_mat <- matrix(NA,2,4)
        var_name <- c("mean","skewness")
        for (i in 1:2){
                ct_tmp <- ct4[ct4$period==i,]
                for (k in 1:2){
                        cor_tmp <- cor.test(log(ct_tmp$rt_est),ct_tmp[,var_name[k]],
                                            use="na.or.complete",
                                            method="spearman")
                        correlation_mat[k,2*(i-1)+1] <- cor_tmp$est
                        correlation_mat[k,2*i] <- cor_tmp$p.value
                }
        }
        
        ### OUTPUT 2: CORRELATION BETWEEN CT-RT ### 
        correlation_mat <- round(correlation_mat,3)
        
        training <- ct4[ct4$date%in%training_period,]
        ### OUTPUT 3: DESCRIPTION OF TRAINING PERIOD ###
        epi_tmp <- epi_plot(training)
        p_training <- epi_tmp[[1]] + 
                annotate(geom="text", x= training$date[round(nrow(training)/3)],
                         y=(epi_tmp[[2]])*0.85,  
                         label=paste0("Duration: ",nrow(training),"d"),
                         fontface="bold",size=3)
        
        ### build model
        lm.used <- lm(log(rt_est)~mean+skewness_imputed,data=training)
        
        rt_compare <- ct4
        est <- exp(predict(lm.used,rt_compare,interval = "prediction"))
        rt_combined <- cbind(rt_compare,est) %>% 
                filter(!is.na(fit))# re-estimated Rt and prediction interval
        
        ### OUTPUT 4.1: CONSISTENCY BETWEEN 2 RT AND SIMULATION TRUTH###
        # i. Spearman rank correlation coefficients (rho)
        # ii. directional consistency
        ## add simulation truth
        rt_truth <- sim_dyn$seir_outputs
        rt_with_real <- merge(rt_combined,rt_truth[,c("step","Rt")],
                              by.x = "sampled_time",by.y = "step",all.x = T) %>%
                filter(period != 0 & !is.na(Rt)) %>%
                mutate(real_cut = factor(cut(Rt,c(0, 0.5, 1, 1.5, 10))))
        
        ### OUTPUT 5: COMPARE WITH "SIMULATION TRUTH" ###
        rt_real_plot <- rt_with_real %>%
                dplyr::select(date,rt_est,fit) %>%
                reshape2::melt(.,id.vars = "date") %>%
                merge(.,rt_with_real[,c("date","real_cut")]) 
        levels(rt_real_plot$variable) <- c("Incidence-based","Ct-based")
        
        #
        p_out2_tmp <- rt_real_plot %>% 
                filter(!is.na(real_cut)) %>%
                ggplot() + geom_boxplot(aes(x = real_cut,y = value,
                                            fill = variable),width = 0.3) + 
                geom_hline(yintercept = 1,
                           linetype = 'dashed',
                           size = 1,
                           color = 'grey') +
                scale_y_continuous(name = 'Estimated Rt',
                                   limits = c(0, 3),
                                   breaks = seq(0, 3, 1),
                                   expand = c(0, 0)) +
                scale_x_discrete(name = 'Simulated truth',
                                 expand = c(0.01, 0.01),
                                 labels = c('<0.5', '0.5-1.0', '1.0-1.5', '>1.5')) +
                scale_fill_manual(name = " ",
                                  values = c('#1b85b8', '#c3cb71')) 
        if (add_legend == T){
                p_out2 <- p_out2_tmp + 
                        theme(legend.position = c(0.75, 0.95)) 
        } else {
                p_out2 <- p_out2_tmp + 
                        theme(legend.position = "none")
        }
        
        p_out3 <- rt_plot(rt_with_real,time = "all", add_epi = T)
        
        #### OUTPUT #####
        return(list(
                # Epi-curve and smoothing Ct (grid.draw() to call the plot)
                trend=p, 
                # correlation between Ct and Rt (training vs. testing)
                correlation=correlation_mat,
                plot0=p_training,
                # model parameters
                model_out=result[[2]],
                # comparison with simulation truth
                plot2=p_out2,
                # comparison between 2 Rt,s
                plot3=p_out3,
                # output dataframe with all Rts
                df_out=rt_with_real)) 
}

###
#####

## end of script

#####