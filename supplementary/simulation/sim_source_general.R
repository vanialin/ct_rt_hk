
## load packages
require(devtools)
#devtools::install_github("jameshay218/lazymcmc")
#devtools::install_github("jameshay218/virosolver")

library(lazymcmc)
library(virosolver)
library(odin)
library(dde)
library(ggplot2)
library(EpiNow2)
library(dplyr)
library(tidyr)
library(reshape2)
library(reprex)
library(tidyverse)
library(patchwork)
library(extraDistr)
library(ggthemes)
library(data.table)
library(mgcv)
library(e1071)
library(extraDistr)
library(MASS)
library(EnvStats)
require(lubridate)
require(gridExtra)
require(chron)
library(purrr)

path_linelist <- paste0(path,"linelists/")
path_rt <- paste0(path_linelist,"output_rt/")
path_plot <- paste0(path,"plots/")

###
source(paste0(path,"sim_funcs_others.R"))
source(paste0(path,"sim_funcs_all.R"))

#
#

#### read in par
model_pars <- read.csv(paste0(path,"partab_seir_switch_model_hk.csv"))
pars <- model_pars$values
names(pars) <- model_pars$names

## Simulated population and time period
population_n <- 7500000
times <- 0:200
## Extend to account for delays
times_extended <-c(times,max(times):(max(times)+50)) 

##
#### end