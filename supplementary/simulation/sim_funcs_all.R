#####
## all sourced functions
#####

library(purrr)
library(EnvStats)

###
#' 1. simulate observed cases 
#'adapted from Hay
#'updated distribution of incubation period and reporting delay

simulate_infected_cases <- function(
        incidence, times, symp_frac=0.6,
        population_n,
        incu_period_par1=log(5.2),incu_period_par2=log(3.9),
        conf_delay_par1=1.83,conf_delay_par2=0.43){
        
        ## Number of individuals who are not infected in the population
        not_infected <- population_n-sum(incidence)
        ## Create a tibble with infection times and incidence
        inc_dat <- tibble(infection_time=c(NA,times),inc=c(not_infected,incidence))
        ## Duplicates rows in inc_dat according to inc
        inc_dat <- inc_dat %>% uncount(inc)
        
        inc_dat <- inc_dat %>%
                mutate(i=1:n()) %>%
                ## If infection time is missing, assign 0 to is_infected, otherwise assign 1
                mutate(is_infected = ifelse(is.na(infection_time), 0, 1)) %>%
                ## Is this individual going to be symptomatic? Uses a binomial distribution
                mutate(is_symp=ifelse(is_infected, rbinom(n(), 1, symp_frac), 0)) %>%
                ## Symptom onset time from a log normal distribution
                mutate(incu_period=ifelse(is_infected & is_symp, rlnormTrunc(n(), 
                                                                             incu_period_par1, 
                                                                             incu_period_par2,
                                                                             min=0,
                                                                             max=25), NA),
                       onset_time=infection_time+floor(incu_period)) %>%
                ## Confirmation time from a gamma distribution (fitted using HK data)
                mutate(confirmation_delay=floor(rgamma(n(), 
                                                       shape = conf_delay_par1, 
                                                       rate = conf_delay_par2)))
        inc_dat
}

###
#' 2. simulate viral trajectory for each infected case
#'adapted from Quilty
#'following assigned incubation period and reporting delay from already simulated values

approx_sd <- function(x1, x2){
        (x2-x1) / (qnorm(0.95) - qnorm(0.05) )
}

get_indiv_trajectory <- function(df){
        #simulate CT trajectories
        set.seed(1)
        traj <- df %>% 
                mutate(u = runif(n(),0,1)) %>%
                # duration from Cevik et al. 2020 Lancet Microbe
                mutate(start=0,
                       end=qnormTrunc(p = u, mean=17, sd=approx_sd(15.5,18.6), min = 0),
                       # we peak -3 ~ 0 d before illness onset (from Jones 2021 Science; He 2020 Nat Med)
                       # +runif(n(),-3,0)
                       onset_t=incu_period) %>% 
                pivot_longer(cols = c(start,end,onset_t),
                             values_to = "x") %>% 
                # peak CT taken from Kissler et al. 2021 PLOS Biology
                # also consistent with Jones 2021 Science
                # which is also quite consistently with the observed highest value in HK data
                mutate(y=case_when(name=="start"   ~ 40,
                                   name=="end"     ~ 40,
                                   name=="onset_t" ~ rnorm(n=n(),mean=22.3,sd=4.2))) %>%  
                ### we get individual trajectories here
                nest(data = c(name,x,y)) %>%  
                dplyr::mutate(
                        # Perform loess calculation on each individual 
                        m  = purrr::map(data, ~splinefunH(x = .x$x, y = .x$y,
                                                          m = c(0,0,0))),
                        pred = purrr::map(.x = m, 
                                          ~data.frame(x = seq(0, 40, by = 1)) %>%
                                                  mutate(y = .x(x)))) %>%
                unnest(pred) %>% 
                rename(infect_to_test=x,ct_value=y) %>%
                dplyr::select(-c(u,data,m)) 
        ##
        ## get the Ct value at sampling
        df_out <- df %>%  
                # here is the actual delay between infection and testing
                mutate(infect_to_test=(floor(incu_period)+confirmation_delay)) %>%
                left_join(traj[,c('i','infect_to_test','ct_value')],by = c("i","infect_to_test")) %>%
                mutate(ct_value=ifelse(is.na(ct_value),40,ct_value))
        
        return(list(traj,df_out))
}

##################
##################

###
#' 3. to get incidence-based Rt from EpiNow2
#' can simply input the observed cases; deconvolution done by the function itself
get_inc_Rt <- function(ct_list0,n1,n2){
        start_time <- Sys.time()
        start_date <- as.Date("2020-07-01") # day 0 for simulation
        # create df for case-counting
        reported_cases <- ct_list0 %>% group_by(confirmed_time) %>%
                # in simulation, confirmation is on the same date as sampling
                summarize(confirm=n()) %>% rename(date=confirmed_time) 
        reported_cases$date <- start_date+reported_cases$date
        
        # estimate Rt
        out <- epinow(reported_cases = reported_cases, 
                      generation_time = generation_time,
                      delays = delay_opts(incubation_period, reporting_delay),
                      rt = rt_opts(prior = list(mean = 2, sd = 2)),
                      gp = gp_opts(basis_prop = 0.2),
                      CrIs = 0.95,
                      stan = stan_opts(cores=4,control=list(adapt_delta=0.95)),
                      horizon = 14, 
                      target_folder = "results",
                      logs = file.path("logs", Sys.Date()),
                      return_output = TRUE, 
                      verbose = TRUE)
        
        # 4. check Rt plots and extracted numbers
        rt <- summary(out, type = "parameters", params = "R")
        write.csv(rt,paste0(path_rt,"rt_obs",n1,"_scenario",n2,".csv"),row.names = F)
        print(paste0("Total time cost = ",
                     round((Sys.time() - start_time),2)," h."))
}

###
#' 4.1. merging Ct and Rt dataframe (for result function)
#' daily mean value and skewness
merge_Ct_Rt <- function(rt0,ct0){
        ct1 <- ct0 %>% filter(ct_value < 40) %>% 
                group_by(sampled_time) %>% 
                summarise(count=n(),
                          mean=mean(ct_value),
                          skewness=skewness(ct_value)) %>%
                ungroup() %>%
                right_join(ct0%>%group_by(infection_time) %>% 
                                   summarise(infect_count=n()) %>% ungroup() %>% 
                                   rename(sampled_time=infection_time))
        
        complete_days <- seq(0,max(ct1$sampled_time),1)
        dummy_row <- tibble(sampled_time=complete_days[!complete_days%in%ct1$sampled_time],
                            count=0,mean=NA,skewness=NA,infect_count=0)
        ct2 <- bind_rows(ct1,dummy_row) %>% arrange(sampled_time)
        
        start_date <- as.Date("2020-07-01")
        ct3 <- ct2 %>% 
                mutate( ## impute skewness
                        skewness_imputed = sapply(ct2$sampled_time, function(d){
                                mean(ct2$skewness[ct2$sampled_time %in% c((d-7):(d-1))], 
                                     na.rm = TRUE)
                        }),
                        date=start_date+sampled_time
                ) 
        ct3$skewness_imputed[!is.na(ct3$skewness)] <- ct3$skewness[!is.na(ct3$skewness)]
        
        rt1 <- rt0 %>% mutate(date=as.Date(date)) %>%
                dplyr::select(date,mean,lower_95,upper_95) %>%
                rename(rt_est=mean,rt_lb=lower_95,rt_ub=upper_95)
        ct4 <- merge(ct3,rt1,by.x = "date",by.y = "date")
        return(ct4)
}

###
#' 4.2. for selecting best-fit training period (for result function)
#' build on estimated incidence-based Rt and simulated daily Ct
select_training_period <- function(df){
        # input data should be merged data (with Ct and inc-Rt)
        initial_start <- df$date[which(df$count>15)][1]
        interval <- 2:6*10 
        r.sqr.mat <- matrix(NA,40,5)
        for (i in 1:40){ # 40 starts
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

###
#' IMPORTANT
#' 5. this output fitted results ***
evaluate_daily_funcs <- function(rt0,ct0,sim_dyn){
        # sim_dyn = seir_dynamics (sourced, simulation truth)
        start_date <- as.Date("2020-07-01") # day 0 for simulation
        #rt0
        # get summary df
        ct4 <- merge_Ct_Rt(rt0,ct0)
        
        ##
        # define training period
        result <- select_training_period(ct4)
        training_period <- seq(result$period[1],
                               result$period[2],1)
        
        ct4$period <- ifelse(ct4$date%in%training_period,1,
                             # 1-training; 2-testing
                             ifelse(ct4$date<min(training_period),0,2))
        
        training <- ct4[ct4$date%in%training_period,]
        
        ### build model
        lm.used <- lm(log(rt_est)~mean+skewness_imputed,data=training)
        
        rt_compare <- ct4
        est <- exp(predict(lm.used,rt_compare,interval = "prediction"))
        rt_combined <- cbind(rt_compare,est) %>% 
                filter(!is.na(fit)) # re-estimated Rt and prediction interval
        
        ## add simulation truth
        rt_truth <- sim_dyn$seir_outputs
        rt_with_real <- merge(rt_combined,rt_truth[,c("step","Rt")] %>% 
                                      mutate(date=start_date+step),
                              by.x = "date",by.y = "date",all.x = T)  %>%
                filter(period!=0)
        
        #### OUTPUT #####
        return(rt_with_real)
}

##################
##################

### 
#' 6. for plotting results
#' (1) epi-curve (by reporting) and Rt curves
fit_plot <- function(df, period,input_range,add_legend,panel){
        # input data is the dataframe generated from evaluate_daily_funcs()
        # with both incidence-based and Ct-based Rt
        y.max = 4
        
        if(period == 'training') {
                data_used = df %>%
                        filter(period == 1)
                dates = range(data_used$date)
                input_range <- NULL
        }
        
        if(period == 'testing') {
                data_used = df %>%
                        filter(date>as.Date(input_range[1])&
                                       date<as.Date(input_range[2]))
                        #filter(period == 2)
                dates = range(data_used$date)
                add_legend = F
                panel = NULL
        }
        
        ### as the main purpose is to compare the consistency between estimates
        ### we did not truncate anything in this function
        ##
        
        y_upper <- round(max(df$count)+150,digits = -2)
        times_to_multiply <- y_upper/4
        
        p = ggplot(data = data_used %>%
                           mutate(upr = ifelse(upr > y.max, y.max, upr))) +
                
                #epi curve 
                geom_bar(data=data_used %>% filter(count > 30),
                         aes(x = date,
                             y = count),
                         stat = 'identity',
                         fill = '#bfbbbb') +
                geom_bar(data=data_used %>% filter(count <= 30),
                         aes(x = date,
                             y = count),
                         stat = 'identity',
                         fill = '#b6d2d7') +
                
                #case rt
                geom_polygon(data = tibble(
                        date_new = c(data_used$date, 
                                     rev(data_used$date)),
                        y_new = c(data_used$rt_lb,
                                  rev(data_used$rt_ub))) %>%
                                mutate(y_new = ifelse(y_new > y.max, y.max, y_new)),
                        aes(x = date_new,
                            y = y_new*times_to_multiply,
                            fill = 'empirical'),
                        color = NA,
                        alpha = 0.2) +
                
                geom_line(data = data_used ,
                          aes(x = date,
                              y = rt_est*times_to_multiply,
                              color = 'empirical'),
                          size = 1) +
                #truth
                geom_line(aes(x = date,
                              y = Rt*times_to_multiply,
                              color = 'truth'),
                          size = 1) +
                #ct rt
                geom_point(aes(x = date,
                               y = fit*times_to_multiply,
                               color = 'predicted'),
                           size = 2) +
                
                geom_segment(aes(x = date,
                                 y = lwr*times_to_multiply,
                                 xend = date,
                                 yend = upr*times_to_multiply,
                                 color = 'predicted'),
                             size = 0.8) +
                
                scale_y_continuous(name = 'Cases',
                                   limits = c(0, y_upper),
                                   expand = c(0, 0),
                                   breaks = seq(0, y_upper, length.out = 6),
                                   position = 'right',
                                   sec.axis = sec_axis(~./times_to_multiply, 
                                                       name = 'Rt')) +
                
                scale_x_date(name = 'Date',
                             limits = c(dates[1]-1,dates[2]+1),
                             date_breaks = "10 days", 
                             date_labels = "%d/%m/%y",
                             expand = c(0.01, 0.01))  +
                
                geom_hline(yintercept = 1*times_to_multiply,
                           linetype = 'dashed',
                           size = 1,
                           color = 'black') 
        if (add_legend == T) {
                p <- p + theme(legend.position = c(0.3,0.95),
                              legend.title = element_text(size = 16, face = 'bold')) 
        } else {
                p <- p + theme(legend.position = "none")
        }
        
        if (!is.null(panel)){
                p <- p + labs(title = panel) +
                        theme(plot.title = element_text(size=20, face = 'bold'))
        } else {
                p <- p + labs(title = "") 
        }
        
        if(period == 'training'){
                
                p = p +
                        scale_fill_manual(name = NULL,
                                          values = c(empirical = '#d9534f',
                                                     predicted = '#428bca',
                                                     truth = "#400000") ,
                                          labels = c('Incidence-based Rt', 
                                                     'Ct-based Rt',
                                                     'Simulation truth'), 
                                          guide = "none") +
                        scale_color_manual(name = NULL,
                                           values = c(empirical = '#d9534f',
                                                      predicted = '#428bca',
                                                      truth = "#400000") ,
                                           labels = c('Incidence-based Rt', 
                                                      'Ct-based Rt',
                                                      'Simulation truth'))
                
        } else {
                
                p = p +
                        scale_fill_manual(name = NULL,
                                          values = c(empirical = '#d9534f',
                                                     predicted = '#428bca',
                                                     truth = "#400000"),
                                          labels = c('Incidence-based Rt', 
                                                     'Ct-based Rt',"Simulation truth"), 
                                          guide = "none") +
                        scale_color_manual(name = NULL,
                                           values = c(empirical = '#d9534f',
                                                      predicted = '#428bca',
                                                      truth = "#400000"),
                                           labels = c('Incidence-based Rt', 
                                                      'Ct-based Rt',"Simulation truth"),
                                           guide="none")
                
                
                
        } 
        
        return(p)
        
}
#
#
#' (2) boxplot (compared with simulation truth)
boxplot_func <- function(df,add_legend,rho_ct){
        df1 <- df %>% mutate(real_cut = factor(cut(Rt,c(0, 0.5, 1, 1.5, 10))))
        
        if (add_legend == T){
                legend_pos <- c(0.75, 0.95)
        } else {
                legend_pos <- "none"
        }
        
        ### COMPARE WITH "SIMULATION TRUTH" ###
        df_for_plot <- df1 %>%
                dplyr::select(date,rt_est,fit) %>%
                reshape2::melt(.,id.vars = "date") %>%
                merge(.,df1[,c("date","real_cut")]) 
        levels(df_for_plot$variable) <- c("Incidence-based","Ct-based")
        
        #
        p_out2_tmp <- df_for_plot %>% 
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
                scale_x_discrete(name = 'Simulation truth',
                                 expand = c(0.01, 0.01),
                                 labels = c('<0.5', '0.5-1.0', '1.0-1.5', '>1.5')) +
                scale_fill_manual(name = " ",
                                  values = c('#d9534f','#428bca')) +
                theme(legend.position = legend_pos) +
                annotate("text", x = 1.6, y = 2.8, 
                         label = paste("rho[ct] == ", rho_ct), parse = TRUE,size=4) +
                labs(title = "") 
        return(p_out2_tmp)
}
#
#
#' (3) Spearman Rho and directional consistency
#' also as compared with simulation truth
Rt_consistency <- function(df,rt1,rt2){
        # rt1 - incidence-based Rt
        # rt2 - Ct-based Rt
        var_name <- c(rt1,rt2)
        cor.rt <- matrix(NA,2,4)
        for (i in 1:2){
                test_tmp <- cor.test(log(df[,var_name[i]]),log(df$Rt),
                                     method = 'spearman', use="na.or.complete",
                                     conf.level = .95)
                cor.rt[i,1:2] <-
                        c(round(test_tmp$est,2),
                          ifelse(test_tmp$p.value<0.001,"<0.001",
                                 round(test_tmp$p.value,3)))
                cor.rt[i,3] <- 
                        round(sum(abs(log(df[,var_name[i]])-log(df$Rt)))/nrow(df),3)
                cor.rt[i,4] <-  
                       paste0(round((nrow(df[((df[,var_name[i]]-1)*(df$Rt-1))>0,])/
                                       nrow(df))*100,1),"%")
        }

        colnames(cor.rt) <- c("rho","p-value","mae","prop.")
        return(cor.rt)
}
##
#####

## end of script

#####