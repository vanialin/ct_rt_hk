#------------
# simulation
# simulate all linelists
# November 2021
#------------

path <- "/Users/vanialam/OneDrive - connect.hku.hk/vanialam/research_vania/epi_wave_2021/program/2021_09_R1/publish (EDIT HERE)/2022_01_R2/"
setwd(path)
source(paste0(path,"sim_source_general.R"))

#########
#########
# simulate "truth" using stochastic SEIR model
# two consecutive waves
# simulate the whole population 
# assign onset, incubation periods and reporting delay

### simulate incidence ----
set.seed(1) 
seir_dynamics <- simulate_seir_wrapper(population_n=population_n,solve_times=times,
                                       pars=pars, switch_model = T,beta_smooth = 0.8)

# simple checking
plot(seir_dynamics$seir_outputs$inc,las=1)
abline(h=1000,col="grey",lty=3)
par(new=T)
plot(seir_dynamics$seir_outputs$Rt,type="l",col="red",axes=F,ylim=c(0,4))
abline(h=1,lty=2,col="red")
abline(h=3,lty=2,col="red")

#save(seir_dynamics,file=paste0(path_linelist,"SEIR_dynamics.Rda"))

### simulate the whole population ----
## if is_infect = 1, assign onset, incubation period and confirmation delay
set.seed(1)
complete_linelist <- simulate_infected_cases(seir_dynamics$incidence,
                                             times=times,population_n=population_n)
write_csv(x=complete_linelist,path=paste0(path_linelist,"complete_linelist.csv"))

symp_linelist <- complete_linelist %>% filter(is_infected==1&is_symp==1)
#write_csv(x=symp_linelist,path=paste0(path_linelist,"symp_linelist.csv"))

#
#

#########
#########
# simulate viral load trajectories for all symptomatic infected individuals
#a <- Sys.time()
set.seed(1)
vl_list <- get_indiv_trajectory(symp_linelist)
save(vl_list,file=paste0(path_linelist,"vl_all_linelist.Rda")) ## this can be large
#Sys.time()-a
vl_traj <- vl_list[[1]] 
### randomly check the simulated trajectories
vl_for_check <- vl_traj %>% filter(i %in% sample(unique(i),200,F))
plot(NA,xlim=c(0,40),ylim=rev(c(10,40)),las=1)
id_check <- unique(vl_for_check$i)
for (n in 1:200){
        lines(vl_for_check$infect_to_test[vl_for_check$i==id_check[n]],
              vl_for_check$ct_value[vl_for_check$i==id_check[n]],col=alpha("grey",.7))
} # checked

### output sampled Ct
vl_full <- vl_list[[2]]
summary(vl_full$ct_value);hist(vl_full$ct_value[vl_full$ct_value<40])
write_csv(x=vl_full,path=paste0(path_linelist,"vl_ob_linelist_full.csv"))

#
#

#########
#########
# simulate each linelinst of detected cases 
# 4 various sets of scenarios 
# get detected cases under each scenario

#### scenarios
set.seed(1)  #### apply for ALL scenarios below ****
### scenario 1&2 - flat detection at 25% and 10%
case_flat_limited <- NULL
for (i in 1:2){
        case_flat_limited_tmp <- 
                simulate_reporting(vl_full,
                                   frac_report = 0.25-0.15*(i==2), 
                                   solve_times=times) %>% arrange(infection_time)
        case_flat_limited[[i]] <- case_flat_limited_tmp
        #write_csv(x=case_flat_limited_tmp,path=paste0(path_linelist,"vl_obs_scenario",i,".csv"))
}

### scenario 3 - increasing detection (per the case of definition changes)
prob_varying <- 
        tibble(t=times_extended,prob=logistic_func(times_extended,
                                                   start_prob=0.15,
                                                   end_prob=0.6,
                                                   growth_rate=0.2,
                                                   # switch point at start of second wave
                                                   switch_point=pars["t_switch2"]), 
               ver="varying")

case_varying <- 
        simulate_reporting(vl_full,
                           timevarying_prob = prob_varying, 
                           solve_times=times) %>% arrange(infection_time)

#write_csv(x=case_varying,path=paste0(path_linelist,"vl_obs_scenario3.csv"))

### scenario 4 - under detection at second wave
a <- logistic_func(101:120,start_prob=0.5,end_prob=0.05,growth_rate=0.25,
                   # switch point before start of second wave
                   switch_point=100)

detect_prob1 <- c(a,rep(a[length(a)],7),rev(a))

prob_ud <- tibble(t=times_extended,
                  prob=c(rep(0.25,100),detect_prob1,
                         rep(0.25,length(times_extended)-100-length(detect_prob1))),
                  ver="ud")
#
case_ud <-  simulate_reporting(vl_full,timevarying_prob = prob_ud, 
                               solve_times=times) %>% 
        arrange(infection_time)
#write_csv(x=case_ud,path=paste0(path_linelist,"vl_obs_scenario4.csv"))
##
#####

## end of script

#####