#------------
# simulate "truth" using stochastic SEIR model
# two consecutive waves
# simulate the whole population 
# assign onset, incubation periods and reporting delay
#------------

## load packages
require(devtools)
#devtools::install_github("jameshay218/lazymcmc")
#devtools::install_github("jameshay218/virosolver")

library(lazymcmc)
library(virosolver)
library(odin)
library(dde)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(reprex)
library(tidyverse)
library(patchwork)
library(extraDistr)
library(ggthemes)
library(EpiNow2)
library(data.table)
library(mgcv)
library(e1071)
library(extraDistr)
library(MASS)

path <- "/Users/vanialam/OneDrive - connect.hku.hk/vanialam/research_vania/epi_wave_2021/program/2021_09_R1/publish (EDIT HERE)/simulations"
setwd(path)
path_source <- paste0(path,"/sources/")
path_simulate <- paste0(path,"/linelists/")
path_complete <- paste0(path_simulate,"complete_linelist/")
path_observe <- paste0(path_simulate,"observe_linelist/")
path_output <- paste0(path_simulate,"output_rt/")
path_plot <- paste0(path,"/plots/")


###
source(paste0(path_source,"odin_funcs.R"))
source(paste0(path_source,"linelist_sim_funcs.R"))

options(mc.cores = 4)

#
#

#### read in par
model_pars <- read.csv(paste0(path,"/sources/partab_seir_switch_model_hk.csv"))
pars <- model_pars$values
names(pars) <- model_pars$names

## Simulated population and time period
population_n <- 1000000
times <- 0:300
## Extend to account for delays
times_extended <-c (times,max(times):(max(times)+50)) 

### simulate linelist for infection (incidence) ----
set.seed(1) 
seir_dynamics <- simulate_seir_wrapper(population_n=population_n,solve_times=times,
                                       pars=pars, ver="odin",switch_model = T,
                                       beta_smooth = 0.8)

# simple checking
plot(seir_dynamics$seir_outputs$inc,las=1)
abline(h=1000,col="grey",lty=3)
par(new=T)
plot(seir_dynamics$seir_outputs$Rt,type="l",col="red",axes=F,ylim=c(0,4))
abline(h=1,lty=2,col="red")
abline(h=3,lty=2,col="red")

save(seir_dynamics,file=paste0(path_simulate,"/SEIR_dynamics.Rda"))

###
# check distribution for reporting delay
set.seed(1)
a <- rdgamma(1000,shape=5,scale=1)
(dist <- fitdistr(a[a>0],"lognormal"))
exp(dist[[1]])

### simulate the whole population ----
## if is_infect = 1, assign onset, incubation period and confirmation delay
set.seed(1)
complete_linelist <- simulate_observations_wrapper(seir_dynamics$incidence,
                                                   times=times,symp_frac=0.6,
                                                   population_n=population_n,
                                                   conf_delay_par1 = 5,
                                                   conf_delay_par2 = 1)
                 
#write_csv(x=complete_linelist, 
 #         path=paste0(path_complete,"complete_linelist_1.csv"))
          
# simulate the otehr linelist for delayed detection (under limited capacity)
complete_linelist2 <- simulate_observations_wrapper(seir_dynamics$incidence,
                                                    times=times,symp_frac=0.6,
                                                    population_n=population_n,
                                                    conf_delay_par1 = 6, ##
                                                    conf_delay_par2 = 1)
        
#write_csv(x=complete_linelist2, 
 #        path=paste0(path_complete,"complete_linelist_1_delay.csv"))
## continue in "sim_2_each_linelist"
#####
#####