### here input functions adapted from source webpage
#' 1. simulate_seir_wrapper
#' 2. logistic_func
#' 3. simulate_reporting
#' 4. odin functions
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

## Subset line list data by detection strategy under symptom-based surveillance 
simulate_reporting <- function(individuals,
                               solve_times, 
                               frac_report=1,
                               timevarying_prob=NULL){
        ## The basic form is that a constant fraction of individuals are reported each day
        n_indivs <- nrow(individuals)
        sampled_individuals <- individuals
        
        ## If a flat reporting rate
        if(is.null(timevarying_prob)){
                ## Symptomatic based surveillance. Observe individuals based on symptom onset date
                ## Subset by symptomatic individuals, observe some fraction of these at some delay post symptom onset
                sampled_individuals <- sampled_individuals %>% 
                        filter(is_symp==1) %>% 
                        sample_frac(frac_report) %>%
                        mutate(sampled_time=onset_time+confirmation_delay,
                               ## We sample individuals some number of days after symptom onset
                               confirmed_time=sampled_time)
                ## Reporting rate varies over time
        } else {
                ## If you have symptom onset on day t, there is a probability that you will be observed
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
        
        return(sampled_individuals=sampled_individuals)
        
}

###

#### for simulate_seir_wrapper() below
## Adapted from
## http://epirecip.es/epicookbook/chapters/sir-stochastic-discretestate-discretetime/r_odin
seir_generator <- odin::odin({
        ## Core equations for transitions between compartments:
        update(S) <- S - n_SE
        update(E) <- E + n_SE - n_EI
        update(I) <- I + n_EI - n_IR
        update(R) <- R + n_IR
        update(inc) <- n_SE
        
        ## Individual probabilities of transition:
        p_SE <- 1 - exp(-beta * I / N) # S to E
        p_EI <- 1 - exp(-sigma) # E to I
        p_IR <- 1 - exp(-gamma) # I to R
        
        ## Draws from binomial distributions for numbers changing between
        ## compartments:
        n_SE <- rbinom(S, p_SE)
        n_EI <- rbinom(E, p_EI)
        n_IR <- rbinom(I, p_IR)
        
        ## Total population size
        N <- S + E + I + R
        
        ## Initial states:
        initial(S) <- S_ini
        initial(E) <- E_ini
        initial(I) <- I_ini
        initial(R) <- 0
        initial(inc) <- 0
        
        ## User defined parameters - default in parentheses:
        S_ini <- user(1000)
        E_ini <- user(0)
        I_ini <- user(1)
        beta <- user(0.2)
        sigma <- user(0.15)
        gamma <- user(0.1)
        
}, verbose = FALSE)
#
## Adapted from 
## http://epirecip.es/epicookbook/chapters/sir-stochastic-discretestate-discretetime/r_odin
seir_generator_switch <- odin::odin({
        ## Core equations for transitions between compartments:
        update(S) <- S - n_SE
        update(E) <- E + n_SE - n_EI
        update(I) <- I + n_EI - n_IR
        update(R) <- R + n_IR
        update(inc) <- n_SE
        
        ## Individual probabilities of transition:
        beta <- beta1*(step <= t_switch1) + beta2*(step > t_switch1 && step < t_switch2) + beta3*(step >= t_switch2)
        
        p_SE <- 1 - exp(-beta * I / N) # S to E
        p_EI <- 1 - exp(-sigma) # E to I
        p_IR <- 1 - exp(-gamma) # I to R
        
        ## Draws from binomial distributions for numbers changing between
        ## compartments:
        n_SE <- rbinom(S, p_SE)
        n_EI <- rbinom(E, p_EI)
        n_IR <- rbinom(I, p_IR)
        
        ## Total population size
        N <- S + E + I + R
        
        ## Initial states:
        initial(S) <- S_ini
        initial(E) <- E_ini
        initial(I) <- I_ini
        initial(R) <- 0
        initial(inc) <- 0
        
        ## User defined parameters - default in parentheses:
        S_ini <- user(1000)
        E_ini <- user(0)
        I_ini <- user(1)
        beta1 <- user(0.2)
        beta2 <- user(0.2)
        beta3 <- user(0.2)
        sigma <- user(0.15)
        gamma <- user(0.1)
        t_switch1 <- user(25)
        t_switch2 <- user(50)
        
        
}, verbose = FALSE)
#
## Adapted from 
## http://epirecip.es/epicookbook/chapters/sir-stochastic-discretestate-discretetime/r_odin
seir_generator_interpolate <- odin::odin({
        ## Core equations for transitions between compartments:
        update(S) <- S - n_SE
        update(E) <- E + n_SE - n_EI
        update(I) <- I + n_EI - n_IR
        update(R) <- R + n_IR
        update(inc) <- n_SE
        
        ## Individual probabilities of transition:
        p_SE <- 1 - exp(-beta * I / N) # S to E
        p_EI <- 1 - exp(-sigma) # E to I
        p_IR <- 1 - exp(-gamma) # I to R
        
        ## Draws from binomial distributions for numbers changing between
        ## compartments:
        n_SE <- rbinom(S, p_SE)
        n_EI <- rbinom(E, p_EI)
        n_IR <- rbinom(I, p_IR)
        
        ## Total population size
        N <- S + E + I + R
        
        ## Initial states:
        initial(S) <- S_ini
        initial(E) <- E_ini
        initial(I) <- I_ini
        initial(R) <- 0
        initial(inc) <- 0
        
        ## User defined parameters - default in parentheses:
        beta = interpolate(betat,betay,"constant")
        output(beta) = beta
        betat[] = user()
        betay[] = user()
        dim(betat) <- user()
        dim(betay) <- length(betat)
        
        S_ini <- user(1000)
        E_ini <- user(0)
        I_ini <- user(1)
        sigma <- user(0.15)
        gamma <- user(0.1)
        
        
}, verbose = FALSE)
##
#####

## end of source file

#####