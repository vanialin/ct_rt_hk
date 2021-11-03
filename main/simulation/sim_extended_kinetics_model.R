#-----------
# Pre-script
# fit kinetics model 
# based on expected parameters
# from HK data
#-----------

## load packages
library(tidyverse)
library(data.table)
library(coda)
library(MASS)
library(doParallel)
library(ggpubr)
library(rethinking)
library(extraDistr)
library(lazymcmc)
library(patchwork)
library(lhs)
library(optimx)
library(ggthemes)
library(mgcv)

path <- "/Users/vanialam/OneDrive - connect.hku.hk/vanialam/research_vania/epi_wave_2021/program/2021_09_R1/publish (EDIT HERE)/simulations"
setwd(path)
path_source <- paste0(path,"/sources/")
##
#source(paste0(path_source,"priors.R"))

parTab <- read.csv(paste0(path_source,"partab_for_optim.csv"))
data <- read.csv(paste0(path_source,"data_ct.csv"))

Mode <- function(x) { # get modal value
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
}

##
ggplot(data) + # if change into log-scale the slope diff will be more significant
        geom_point(aes(x=test.to.onset, y=log(ct.value)), alpha=.4) +
        geom_smooth(aes(x=test.to.onset,y=log(ct.value))) +
        geom_vline(xintercept = c(12,16),linetype = 'dashed',size = .4,color = 'red') + 
        scale_y_reverse() + coord_cartesian(xlim=c(0,50), ylim=log(c(40,10)))
#
# Modal Ct value at a = teclipse (0) + tpeak (5) + tswitch (12-16 in HK data)
parTab[parTab$names == "level_switch","values"] <- # value = 35
        Mode(data$ct.value[data$test.to.onset%in%12:16])
parTab[parTab$names == "t_switch","value"] <- 13 # from shedding paper?
ct_lower <- range(data$ct.value,na.rm=TRUE)[1]
# modal Ct value at day 30 **post infection**
ct_upper <- round(Mode(data$ct.value[data$test.to.onset%in%24:26])) # value = 33


#
#

##### MODEL -----
## 1) pre-requisites
parTab[parTab$names == "LOD","values"] <- 3

assumed_incu_period <- 5

## Detection probability as a function of days since symptom onset
## Add assumed incubation period to get days since infection
to_fit_probs <- read_csv(paste0(path_source,"borremands_urt_pos.csv"))
ages_observed <- 
        c(-assumed_incu_period, to_fit_probs$time_since_onset) + assumed_incu_period
ages <- ages_observed

desired_probs <- (c(0,to_fit_probs$pos)/100)
desired_probs <- pmax(desired_probs, 0)

targeted_dat <- tibble(age=ages_observed, prob=desired_probs)

pars_hinge2 <- parTab$values
names(pars_hinge2) <- parTab$names

#

## 2) cost function
cost_function_hinge <- function(explore_pars, use_ct_cost=TRUE, ct_cost_weight=1, 
                                ct_lower=ct_lower,ct_upper=35){
        pars1 <- pars_hinge2
        
        ## Pull out the parameters being optimized
        pars1["t_switch"] <- explore_pars[1]
        pars1["obs_sd"] <- explore_pars[2]
        pars1["viral_peak"] <- explore_pars[3]
        pars1["level_switch"] <- explore_pars[4]
        #pars1["wane_rate2"] <- explore_pars[5]
        pars1["prob_detect"] <- explore_pars[6]
        pars1["sd_mod"] <- explore_pars[7]
        
        ## Find Cts at peak and at day 30
        peak <- viral_load_func(pars1,pars1["desired_mode"]+pars1["tshift"],convert_vl=FALSE)
        day30 <- viral_load_func(pars1,30,convert_vl=FALSE)
        
        ## Calculate cost function for peak and day 30 Cts
        cost1 <- 0
        cost3 <- 0
        ## Also include highest recorded Ct
        ## Lower 1% of distribution should be ct_upper at peak and
        ## upper 1% should be ct_lower at day 30
        if(use_ct_cost){
                cost1 <- (ct_lower - qgumbel(0.01,peak,pars1["obs_sd"]))^2
                cost3 <- (ct_upper - qgumbel(0.01,day30,pars1["obs_sd"]*pars1["sd_mod"]))^2
        }
        ## Calculate proportion detectable curve
        ## Get proportion detectable over time and fit to observations
        vls <- viral_load_func(pars1, 1:max(ages), convert_vl=FALSE)
        obs <- c(0,prop_detectable(ages[ages > 0], pars1, vls))
        costs2 <- (obs-desired_probs)^2
        sum(costs2 + cost1*ct_cost_weight + cost3)
}

#

## 3) OPTIM
set.seed(1)
fit_hinge2 <- optimx(c(10,5.5,30,40,40,0.2,0.8), cost_function_hinge,
                     lower=c(0,5,0,30,15,.000001,0.25),
                     upper=c(25,10,40,40,100,1,1),
                     method="L-BFGS-B",
                     use_ct_cost=TRUE, ct_cost_weight=1,
                     ct_lower=ct_lower,ct_upper=ct_upper)

pars1 <- pars_hinge2

pars1["t_switch"] <- as.numeric(fit_hinge2[1,1])
pars1["obs_sd"] <- as.numeric(fit_hinge2[1,2])
pars1["viral_peak"] <- as.numeric(fit_hinge2[1,3])
pars1["level_switch"] <- as.numeric(fit_hinge2[1,4])
pars1["prob_detect"] <- as.numeric(fit_hinge2[1,6])
pars1["sd_mod"] <- as.numeric(fit_hinge2[1,7])
pars1["t_unit"] <- 1

#

# write into par for modelling
parTab$values <- pars1
write_csv(parTab,paste0(path_source,"partab_fitted_hk.csv"))

##
# updated input pars
orig_par <- read.csv(paste0(path_source,"partab_seir_switch_model_v2.csv"),as.is = T)
change_par <- orig_par
for (i in 1:nrow(change_par)){
        if (change_par$names[i]%in%parTab$names){
                change_par$values[i] <- 
                        parTab$values[parTab$names==change_par$names[i]]
        }
}
write_csv(change_par,paste0(path_source,"partab_seir_switch_model_hk.csv"))
##
#####

## end of script

#####