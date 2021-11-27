### here put those directly input from others' codes
## Functions here:
#' 2. simulate_seir_wrapper
#' 3. logistic_func
#' 4. simulate_reporting
#' adapted to allow sharper epi waves 
#' on November 17, 2021

## Wrapper for SEIR model simulation, deterministic or stochastic
## INPUTS: 
##      1. population_n: size of population to simulate
##      2. solve_times: vector of times to solve model over
##      3. pars: named vector of SEIR model parameters
##      4. version: either "ode" for determinsitic version or "odin" for stochastic
##      5. switch_model: if TRUE, uses the SEIR model with 2 switch points in transmission intensity
##      6. beta_smooth: spar parameter fed to smooth.spline, used to smooth the time-varying beta
## OUTPUTS: 
##      1. Data frame of the SEIR solution
##      2. Vector of absolute incidence per time point
##      3. Overall probability of infection
##      4. Plot of incidence over time
##      5. Average growth rates over different time periods
##      6. Plot of average growth rates over different time periods
simulate_seir_wrapper <- function(population_n, solve_times, pars, version="ode",switch_model=FALSE,beta_smooth=0.8){
        ## Choose which version of the model to solve
        if(version == "ode"){
                ## Deterministic SEIR model
                epidemic_process <- simulate_seir_process(pars,solve_times,N=population_n,switch_model=switch_model)
                
                ## Widen solution and extract key variables
                res <- epidemic_process$seir_outputs %>% pivot_wider(names_from="variable",values_from="value")
                res <- res %>% rename(step=time)
                incidence <- epidemic_process$per_cap_incidence
                overall_prob <- epidemic_process$overall_prob_infection
        } else {
                ####################################################
                ## Stochastic model
                ####################################################
                if(switch_model){
                        gamma1 <- 1/pars["infectious"]
                        sigma1 <- 1/pars["incubation"]
                        beta1 <- pars["R0_1"]*gamma1
                        beta2 <- pars["R0_2"]*gamma1
                        beta3 <- pars["R0_3"]*gamma1
                        beta4 <- pars["R0_4"]*gamma1
                        I0 <- ceiling(pars["I0"]*population_n)
                        ## Odin stochastic SEIR model generator
                        #seir <- seir_generator_switch(beta1=beta1,beta2=beta2,beta3=beta3,
                        #                              sigma=sigma1,gamma=gamma1,
                        #                              S_ini=population_n-I0,I_ini=I0,
                        #                              t_switch1=pars["t_switch1"],t_switch2=pars["t_switch2"])
                        betas <- rep(beta4, length(solve_times))
                        betas[which(solve_times < pars["t_switch3"])] <- beta3
                        betas[which(solve_times < pars["t_switch2"])] <- beta2
                        betas[which(solve_times < pars["t_switch1"])] <- beta1
                        betas <- smooth.spline(betas,spar=beta_smooth)$y
                        seir <- seir_generator_interpolate$new(betat=solve_times,betay=betas,sigma=sigma1,gamma=gamma1,S_ini=population_n-I0,I_ini=I0)
                        
                } else {
                        gamma1 <- 1/pars["infectious"]
                        sigma1 <- 1/pars["incubation"]
                        beta1 <- pars["R0_1"]*gamma1
                        I0 <- ceiling(pars["I0"]*population_n)
                        ## Odin stochastic SEIR model generator
                        seir <- seir_generator$new(beta=beta1,sigma=sigma1,gamma=gamma1,S_ini=population_n-I0,I_ini=I0)
                }
                ## Solve model
                res <- seir$run(solve_times)
                ## Make sure we get a simulation with an outbreak - keep trying until it takes off
                while(max(res[,"I"]) <= I0) res <- seir$run(solve_times)
                res <- as.data.frame(res)
                
                ## Shift for start
                res$step <- res$step + floor(pars["t0"])
                
                ## Dummy rows from pre-seeding
                if(pars["t0"] > 0){
                        dummy_row <- data.frame("step"=0:(floor(unname(pars["t0"]))-1),"S"=population_n,"E"=0,"I"=0,"R"=0,"inc"=0)
                        res <- bind_rows(dummy_row, res)
                }
                res <- res[res$step %in% times,]
                
                ## Get raw incidence and overall probability of infection
                incidence <- res$inc/population_n
                overall_prob <- max(res$R)/population_n
                
                if(!switch_model){
                        res$beta <- pars["R0_1"]/pars["infectious"]
                }
                
                res$Rt <- (res$S/population_n) * res$beta * pars["infectious"]
        } 
        ## Get absolute incidence
        incidence <- incidence * population_n
        
        ## Reshape solution
        res_melted <- res %>% pivot_longer(-step)
        
        ## Compartment plot
        p_compartments <- res_melted %>% 
                filter(name %in% c("S","E","I","R","cumulative_incidence")) %>%
                ggplot() + geom_line(aes(x=step,y=value,col=name)) +
                ylab("Per capita") +
                xlab("Date") +
                theme_bw() +
                theme(legend.position="top")
        
        ## Incidence plot
        p_inc <- ggplot(data.frame(x=solve_times,y=incidence)) +
                geom_line(aes(x=x,y=y),col="red") +
                ylab("True incidence") +
                xlab("Date") +
                theme_bw()
        
        ## Rt plot
        p_rt <- res_melted %>% filter(name == "Rt") %>%
                ggplot() +
                geom_line(aes(x=step,y=value),col="blue") +
                scale_y_continuous(limits=c(0,pars["R0"]+1)) +
                geom_hline(yintercept=1,linetype="dashed") +
                ylab("Rt") +
                xlab("Date") +
                theme_bw()
        
        ## Combine plots
        inc_plot <- p_compartments / p_inc / p_rt
        
        ## Get growth rates
        GR_daily <- log(incidence[2:length(incidence)]/incidence[1:(length(incidence)-1)])
        GR_daily <- ifelse(is.finite(GR_daily), GR_daily, NA)
        GR_daily_dat <- data.frame(t=solve_times[2:length(solve_times)],GR=GR_daily,ver="daily")
        
        ## Get average growth rate over different size windows
        lastdays <- seq(10,50,by=10)
        
        GR_all <- GR_daily_dat
        for(t in seq_along(lastdays)){
                GR_full <- NULL
                lastday <- lastdays[t]
                for (i in (lastday+1):length(solve_times)) {
                        end_index <- i-1
                        start_index <- max(1, (i-lastday))
                        GR_full <- c(GR_full,mean(GR_daily[start_index:end_index], na.rm=TRUE))
                }
                GR_full_dat <- data.frame(t=(lastday+1):length(solve_times), GR=GR_full,ver=as.character(lastday))
                GR_all <- bind_rows(GR_all, GR_full_dat)
        }
        
        ## Get daily growth rate around the peak
        gr_crossover <- GR_all %>% filter(ver == "daily") %>%
                filter(t < 250 & t > 100) %>%
                mutate(abs_gr = abs(GR)) %>%
                filter(abs_gr == min(abs_gr, na.rm=TRUE)) %>% pull(t)
        
        ## Average growth rates
        p_gr <- ggplot(GR_all %>% filter(ver != "daily")) +
                geom_line(aes(x=t,y=GR,col=ver)) +
                geom_hline(yintercept=0,linetype="dashed") +
                geom_vline(xintercept=gr_crossover,linetype="dotted")+
                coord_cartesian(ylim=c(-0.2,0.2)) +
                ylab("Growth rate") +
                xlab("Date") +    
                ggtitle("Average growth rate over different windows") +
                theme_bw()
        
        ## Daily growth rate
        p_gr1 <- ggplot(GR_all %>% filter(ver == "daily")) +
                geom_line(aes(x=t,y=GR),col="black") +
                geom_hline(yintercept=0,linetype="dashed") +
                geom_vline(xintercept=gr_crossover,linetype="dotted")+
                coord_cartesian(ylim=c(-0.2,0.2)) +
                ylab("Growth rate") +
                ggtitle("Daily growth rate") +
                xlab("Date") +    
                theme_bw()
        
        
        list(seir_outputs=res, 
             incidence=incidence,
             overall_prob=overall_prob,
             plot=inc_plot, 
             growth_rates=GR_all,
             growth_rate_p =p_gr1/p_gr)
}


## Test an varying fraction of people per day
## INPUTS: 
##      1. t: times to solve over
##      2. start_prob: y value at start
##      3. end_prob: y value at end
##      4. growth_rate: growth rate of logistic function
##      5. switch_point: x value for the inflection point
## OUTPUTS: 
##      1. a vector giving y values corresponding to t
logistic_func <- function(t,start_prob,end_prob, growth_rate, switch_point=100){
        start_prob + (end_prob-start_prob)/(1 + exp(-growth_rate*(t-switch_point)))
}

## Subset line list data by testing strategy. Options:
#' 1. Sample a random fraction of the population if the only argument is frac_report
#' 2. Sample some random fraction of the population at a subset of time points, specified by timevarying_prob
#' 3. Observe symptomatic individuals with some fixed probability, frac_report if symptomatic is TRUE
#' 4. Observe symptomatic individuals with some time-varying probability, timevarying_prob, if symptomatic is TRUE
#' INPUTS: 
#'      1. individuals: the full line list from the simulation, returned by virosolver::simulate_observations_wrapper
#'      2. solve_times: vector of times at which individuals can be reported
#'      3. frac_report: the overall fraction/probability of individuals who are reported
#'      4. timevarying_prob: a tibble with variables t and prob. This gives the probability of being reported on day t
#'      5. symptomatic: if TRUE, then individuals are reported after developing symptoms. If FALSE, then we take a random cross-section 
#' OUTPUTS: 
#'      1. A tibble with line list data for individuals who were observed
#'      2. A plot of incidence for both observed individuals and the entire simulated population
#'      3. Plot growth rate of cases/infections in the entire population and observed population
simulate_reporting <- function(individuals,
                               solve_times, 
                               frac_report=1,
                               timevarying_prob=NULL, 
                               symptomatic=FALSE){
        ## The basic form is that a constant fraction of individuals are reported each day
        n_indivs <- nrow(individuals)
        sampled_individuals <- individuals
        
        ## If a flat reporting rate
        if(is.null(timevarying_prob)){
                ## Base case, just sampled individuals at random
                if(!symptomatic){
                        ## Sample a fraction of the overall line list and assign each individual a sample time (at random)
                        ## and a confirmation time (with confirmation delay)
                        sampled_individuals <- sampled_individuals %>% 
                                sample_frac(frac_report) %>%
                                group_by(i) %>%
                                mutate(sampled_time=sample(solve_times,n()),
                                       ## We sample but then have to wait for a result
                                       confirmed_time=sampled_time+confirmation_delay) %>%
                                ungroup()
                } else {
                        ## Symptomatic based surveillance. Observe individuals based on symptom onset date
                        ## Subset by symptomatic individuals, observe some fraction of these at some delay post symptom onset
                        sampled_individuals <- sampled_individuals %>% 
                                filter(is_symp==1) %>% 
                                sample_frac(frac_report) %>%
                                mutate(sampled_time=onset_time+confirmation_delay,
                                       ## We sample individuals some number of days after symptom onset
                                       confirmed_time=sampled_time)
                }
                ## Reporting rate varies over time
        } else {
                ## Sample entire population
                if(!symptomatic){
                        indivs <- unique(sampled_individuals$i)
                        indivs_test <- indivs
                        tmp_sampled_indivs <- NULL
                        for(index in 1:nrow(timevarying_prob)) {
                                sample_n <- round(timevarying_prob$prob[index]*length(indivs))
                                sampled_time1 <- timevarying_prob$t[index]
                                sampled_indivs <- sample(indivs_test, sample_n,replace=FALSE)
                                tmp_sampled_indivs[[index]] <- sampled_individuals %>% 
                                        filter(i %in% sampled_indivs) %>% 
                                        mutate(sampled_time=sampled_time1,
                                               confirmed_time = sampled_time + confirmation_delay)
                                indivs_test <- setdiff(indivs_test, sampled_indivs)
                        }
                        sampled_individuals <- do.call("bind_rows", tmp_sampled_indivs)
                        
                        #browser()
                        ## How many individuals are we going to observe by the end?
                        #frac_report_overall <- 1-prod(1-timevarying_prob$prob)
                        #frac_report_overall <- sum(timevarying_prob$prob)
                        ## We will first get the fraction of individuals we will sample over the whole period
                        ## Then, we will change the time-varying reporting probability to the relative number
                        ## sampled on each day
                        #scaled_timevarying_prob <- timevarying_prob$prob/frac_report_overall
                        
                        ## On each time point in timevarying_prob$t, choose some random fraction of the population
                        ## to observe, weighted by scaled_timevarying_prob
                        ## assign a sample time and confirmation time
                        #sampled_individuals <- sampled_individuals %>%
                        #  sample_frac(frac_report_overall) %>%
                        # group_by(i) %>%
                        #  mutate(sampled_time = sample(timevarying_prob$t, n(), replace=TRUE, prob=scaled_timevarying_prob),
                        #        confirmed_time = sampled_time + confirmation_delay) %>%
                        #  ungroup()
                } else {
                        ## This is quite different - if you have symptom onset on day t, there is a probability that you will be observed
                        ## Symptomatic based surveillance. Observe individuals based on symptom onset date
                        sampled_individuals <- sampled_individuals %>% 
                                filter(is_symp==1) %>% ## Subset for symptomatic
                                left_join(timevarying_prob %>% rename(onset_time=t) %>% ## Join with time-varying reporting rate table
                                                  dplyr::select(-ver), by="onset_time") %>%
                                group_by(i) 
                        
                        ## Quicker to vectorize
                        sampled_individuals$is_reported <- rbinom(nrow(sampled_individuals), 1, sampled_individuals$prob) ## Assign individuals as detected or not
                        
                        sampled_individuals <- sampled_individuals %>%
                                filter(is_reported == 1) %>% ## Only take detected individuals
                                mutate(sampled_time=onset_time+confirmation_delay,
                                       ## We sample individuals some number of days after symptom onset
                                       confirmed_time=sampled_time)
                }
        }
        
        ## Plot incidence of infections,  onsets, confirmations and number sampled per day
        ## Get grouped (not line list) subset data
        grouped_dat <- sampled_individuals %>% 
                dplyr::select(i, infection_time, onset_time, sampled_time, confirmed_time) %>%
                pivot_longer(-i) %>% 
                drop_na() %>%
                group_by(name, value) %>% 
                tally() %>%
                rename(var=name,
                       t=value,
                       n=n) %>%
                complete(var, nesting(t),fill = list(n = 0)) %>%
                mutate(ver="Sampled individuals")
        
        ## Grouped entire dataset
        grouped_dat_all <- individuals %>% 
                dplyr::select(i, infection_time, onset_time) %>%
                pivot_longer(-i) %>% 
                drop_na() %>%
                group_by(name, value) %>% 
                tally() %>%
                rename(var=name,
                       t=value,
                       n=n) %>%
                complete(var, nesting(t),fill = list(n = 0)) %>%
                mutate(ver="All individuals")
        
        grouped_dat_combined <- bind_rows(grouped_dat, grouped_dat_all)
        
        ## Plot line list data
        p_all <- grouped_dat_combined %>% ggplot() +
                geom_line(aes(x=t,y=n,col=var)) +
                theme_bw() +
                ylab("Number of individuals") +
                xlab("Time") +
                theme(legend.position="bottom") +
                facet_wrap(~ver,ncol=1, scales="free_y")
        
        ## Get day-by-day growth rate
        grouped_dat_combined <- grouped_dat_combined %>% group_by(var, ver) %>% 
                mutate(gr_daily=log(n/lag(n,1))) %>%
                mutate(gr_daily = ifelse(is.finite(gr_daily), gr_daily, NA)) %>%
                ungroup()
        
        ## Get rolling average growth rate
        growth_rates_all <- expand_grid(grouped_dat_combined, window=seq(10,50,by=10)) %>% 
                arrange(var, ver, window, t) %>%
                group_by(var, ver, window) %>%
                mutate(gr_window=zoo::rollmean(gr_daily, window,align="right",fill=NA)) %>%
                mutate(window = as.factor(window))
        
        p_gr <- ggplot(growth_rates_all) +
                geom_line(aes(x=t,y=gr_window,col=ver)) +
                geom_hline(yintercept=0,linetype="dashed") +
                coord_cartesian(ylim=c(-0.2,0.2)) +
                ylab("Growth rate") +
                xlab("Date") +    
                theme_bw() +
                theme(legend.position="bottom") +
                facet_grid(window~var)
        
        
        list(sampled_individuals=sampled_individuals,
             plot=p_all, 
             plot_gr=p_gr)
}

##
#####

## end of source file

#####