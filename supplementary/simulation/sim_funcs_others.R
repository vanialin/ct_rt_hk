### here input functions adapted from source webpage
#' 1. simulate_seir_wrapper
#' 2. logistic_func
#' 3. simulate_reporting
#' November 17, 2021

## source: https://github.com/jameshay218/virosolver/blob/master/R/simulation_functions.R

## Wrapper for stochastic SEIR model simulation
simulate_seir_wrapper <- function(population_n, solve_times, pars,
                                  switch_model=FALSE,beta_smooth=0.8){
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
                dummy_row <- data.frame("step"=0:(floor(unname(pars["t0"]))-1),
                                        "S"=population_n,"E"=0,"I"=0,"R"=0,"inc"=0)
                res <- bind_rows(dummy_row, res)
        }
        res <- res[res$step %in% times,]
        
        ## Get raw incidence and overall probability of infection
        incidence <- res$inc/population_n

        if(!switch_model){
                res$beta <- pars["R0_1"]/pars["infectious"]
        }
        
        res$Rt <- (res$S/population_n) * res$beta * pars["infectious"]
        ## Get absolute incidence
        incidence <- incidence * population_n
        
        list(seir_outputs=res, 
             incidence=incidence)
}


## Test an varying fraction of people per day
logistic_func <- function(t,start_prob,end_prob, growth_rate, switch_point=100){
        start_prob + (end_prob-start_prob)/(1 + exp(-growth_rate*(t-switch_point)))
}

## Subset line list data by testing strategy. 
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
        
        return(sampled_individuals=sampled_individuals)
        
}

##
#####

## end of source file

#####