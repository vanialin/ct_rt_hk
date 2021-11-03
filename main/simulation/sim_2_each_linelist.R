#------------
# simulate each linelinst of detected cases 
# 9 various sets of scenarios
# assign sampled Ct values for each detected cases
#------------

##
### generate different sets of varying detection probabilities---- 
## to mimic varying detection prob. under symptom-based surveillance
prob_varying_list <- list()
p_start <- c(0.1,0.25,0.5,0.6)
p_end <- c(0.3,1,0.6,0.5)
#
for (i in 1:4){
        prob_varying_list[[i]] <- 
                tibble(t=times_extended,prob=logistic_func(times_extended,
                                                           start_prob=p_start[i],
                                                           end_prob=p_end[i],
                                                           growth_rate=0.1,
                                                           # switch point at start of second wave
                                                           switch_point=160), 
                       ver=paste0("scenario ",i))
}
#
## check probability trend over time
plot(NA,xlim=c(0,350),ylim=c(0.1,1),las=1,xlab="Time",ylab="Prob.")
for (i in 1:4){ # only show the first four scenarios
        # changing from wave 1 = dashed line
        lines(prob_varying_list[[i]]$prob,col=i)
}
abline(v=160,lty=3,col="grey")
abline(v=130,lty=3,col="grey")
abline(v=200,lty=3,col="grey")

#
#

### simulate each linelinst of detected cases ----
## under flat but limited detection or varying detection 
reporting_prob <- c(1,.5,.3,.1)
set.seed(1)
### flat detection ----
for (i in 1:4){
        vl_flat_limited_tmp <- 
                simulate_viral_loads_wrapper(
                        simulate_reporting(complete_linelist,
                                           frac_report = reporting_prob[i],
                                           solve_times=times, 
                                           symptomatic=T)$sampled_individuals %>%
                                arrange(infection_time),
                        pars) 
        write_csv(x=vl_flat_limited_tmp, 
                  path=paste0(path_observe,"vl_obs1_scenario",
                              i,".csv"))
}
#
#
### varying detection ----
for (k in 1:4){
        vl_varying_tmp <- simulate_viral_loads_wrapper(
                simulate_reporting(complete_linelist,
                                   timevarying_prob = 
                                           prob_varying_list[[k]], 
                                   solve_times=times, 
                                   symptomatic=T)$sampled_individual %>% 
                        arrange(infection_time),
                pars)
        write_csv(x=vl_varying_tmp,
                  path=paste0(path_observe,"vl_obs1_scenario",4+k,".csv"))
}
##
## limited detection and therefore also delayed detection  
set.seed(1)
vl_flat_delay_tmp <- 
        simulate_viral_loads_wrapper(
                simulate_reporting(complete_linelist2,
                                   frac_report = 0.3,
                                   solve_times=times, 
                                   symptomatic=T)$sampled_individuals %>%
                        arrange(infection_time),
                pars)            
                        
write_csv(x=vl_flat_delay_tmp,
          path=paste0(path_observe,"vl_obs1_scenario_delay.csv"))

##        
#####

## end of script

#####