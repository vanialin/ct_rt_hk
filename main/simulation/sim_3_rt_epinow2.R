#------------
# estimate Rt (EpiNow2)
# based on cases detected under each scenario
#------------

## load extra packages
require(data.table)
require(EpiNow2)

options(mc.cores = 4)

## ALL RUN in WINDOWS ##

#
#

### delay distributions ----
incubation_period <- 
        get_incubation_period(disease = "SARS-CoV-2", 
                              source = "lauer")
#
# below is the GT and reporting delay dist. for simulated data set (adapted from HAY)
SI_mean <- pars["infectious"] + pars["incubation"]
SI_sd <- sqrt(pars["infectious"]^2 + pars["incubation"]^2)
generation_time <- list(mean=SI_mean,mean_sd=3,sd=SI_sd,sd_sd=3,max=100)
#
set.seed(1)
reporting_delay <- 
        bootstrapped_dist_fit(extraDistr::rdgamma(1000, 5, 1), 
                              max_value = 15,bootstraps = 1)

start_date <- as.Date("2020-07-01") # day 0 for simulation

#
#

### run EpiNow2 for each linelist ----
## each run takes around 10 hours
for (i in 1:9){ 
        start_time <- Sys.time()
        # 1. # input dataframe 
        df <- read.csv(paste0(path_observe,
                              "vl_obs1_scenario",i,".csv")) 
        
        # 2. create df for case-counting
        reported_cases <- df %>% group_by(confirmed_time) %>%
                # in simulation, confirmation is on the same date as sampling
                summarize(confirm=n()) %>% rename(date=confirmed_time) 
        reported_cases$date <- start_date+reported_cases$date
        
        # 3. estimate
        #future::plan("multiprocess") # ALTERNATIVE
        out <- epinow(reported_cases = reported_cases, 
                      generation_time = generation_time,
                      delays = delay_opts(incubation_period, reporting_delay),
                      rt = rt_opts(prior = list(mean = 2.5, sd = 2)),
                      gp = gp_opts(basis_prop = 0.2),
                      CrIs = c(0.5, 0.95),
                      stan = stan_opts(cores=4,control=list(adapt_delta=0.99)),
                      #stan = stan_opts(future = TRUE, max_execution_time = 60 * 30), # ALTERNATIVE
                      horizon = 14, 
                      target_folder = "results",
                      logs = file.path("logs", Sys.Date()),
                      return_output = TRUE, 
                      verbose = TRUE)
        
        # 4. check Rt plots and extracted numbers
        rt <- summary(out, type = "parameters", params = "R")
        write.csv(rt,paste0(path_output,"rt_obs1_scenario",i,".csv"),
                  row.names = F)
        print(paste0("Total time cost = ",
                     round((Sys.time() - start_time),2)," h.")) 
} 
##
#####

## end of script

#####