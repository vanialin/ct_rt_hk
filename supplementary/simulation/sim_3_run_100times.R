#------------
# simulation
# run 100 times 
# get table S6
# November 2021
#------------

path <- "HERE IS WHERE YOU PUT YOUR CODES"
setwd(path)
source(paste0(path,"sim_source_general.R"))

#### components ready ####
### data
vl_full <- read.csv(paste0(path_linelist,"vl_ob_linelist_full.csv"))
load(file=paste0(path_linelist,"SEIR_dynamics.Rda"))
#
### functions
# 1) to sample equal number of cases per day as shown in the representative run
sample_ct_only <- function(sampled_individuals,daily_count){
        indivs_all <- unique(sampled_individuals$i)
        sampled_individuals <- sampled_individuals %>% 
                mutate(sampled_time = onset_time+confirmation_delay,
                       confirmed_time = sampled_time)
        tmp_sampled_indivs <- NULL
        for(d in 1:nrow(daily_count)) {
                indivs <- unique(sampled_individuals$i[sampled_individuals$sampled_time==
                                                               daily_count$sampled_time[d]])
                tmp_sampled_indivs[[d]] <- sampled_individuals %>% filter(i%in%sample(indivs,
                                                                                      daily_count$n[d],
                                                                                      replace=FALSE))
        }
        sampled_individuals_out <- do.call("bind_rows", tmp_sampled_indivs) %>% 
                arrange(infection_time)
}
#

# 2) to generate results as CI
get_ci <- function(x){
        round(c(mean(x),quantile(x,c(.025,.975))),3)
}

check_result <- function(seq,run_time){
        result_all <- matrix(NA,run_time,6)
        train_date <- result_list[[seq]]$date[result_list[[seq]]$period==1] # apply to all
        for (k in 1:run_time){
                set.seed(k)
                ct_renew <- sample_ct_only(vl_full,case_count_list[[seq]])
                df <- merge_Ct_Rt(rt_list[[seq]],ct_renew) %>%
                        right_join(rt_list[[seq]] %>% dplyr::select(date,Rt))
                lm.tmp <- lm(log(rt_est)~mean+skewness_imputed,
                             data=df %>% filter(date%in%train_date))
                est <- exp(predict(lm.tmp,df,interval = "prediction"))
                df1 <- cbind(df,est) %>% 
                        mutate(period=ifelse(df$date%in%train_date,1,
                                             # 1-training; 2-testing
                                             ifelse(df$date<min(train_date),0,2))) %>% 
                        filter(!is.na(fit)) 
                for (j in 1:2){
                        df2 <- df1 %>% filter(period==j)
                        result_all[k,2*(j-1)+1] <- 
                                round(cor.test(log(df2$fit),log(df2$Rt),method = 'spearman',
                                               use="na.or.complete",conf.level = .95)$est,2)
                        result_all[k,2*j] <- 
                                round((nrow(df2[((df2$fit-1)*(df2$Rt-1))>0,])/nrow(df2))*100,1)
                        if (j == 2){
                                df3 <- df2 %>% filter(count>30)
                                result_all[k,5] <- 
                                        cor.test(log(df3$fit),log(df3$Rt),method = 'spearman',
                                                       use="na.or.complete",conf.level = .95)$est
                                result_all[k,6] <- 
                                        (nrow(df3[((df3$fit-1)*(df3$Rt-1))>0,])/nrow(df3))*100
                        }
                }
                
        }
        
        result_out <- apply(result_all, 2, get_ci)
        return(result_out)
}

#
#

#### generate results ####
### read in main line lists
ct_list <- rt_list <- result_list <- 
        case_count_list <- NULL ## the rt_list is a bit different from the one in "sim_2"
for (i in 1:4){
        ct_tmp <- read_csv(paste0(path_linelist,"vl_obs_scenario",i,".csv"))
        rt_tmp <- read_csv(paste0(path_rt,"rt_obs1_scenario",i,".csv"))
        result_list[[i]] <- evaluate_daily_funcs(rt_tmp,ct_tmp,seir_dynamics)
        
        rt_tmp1 <- rt_tmp %>% 
                ## merge with simulation truth in this version
                right_join(seir_dynamics$seir_outputs %>% dplyr::select(step,Rt) %>%
                                   mutate(date=as.Date("2020-07-01")+step))
        
        # this line list is used to sample equal number of cases by date of sampling
        # just as shown in Fig. S10
        case_count_list[[i]] <- ct_tmp %>% group_by(sampled_time) %>%
                summarise(n=n()) %>% ungroup() 
        
        ct_list[[i]] <- ct_tmp
        rt_list[[i]] <- rt_tmp1
}

#
con_result <- NULL
for (i in 1:4){
        con_result[[i]] <- 
                check_result(seq = i,run_time = 100)
}


out <- matrix(NA,4,6)
for (i in 1:4){
        for (k in 1:3){
                out[i,2*(k-1)+1] <- 
                        paste0(round(con_result[[i]][1,2*(k-1)+1],2),
                               "(",round(con_result[[i]][2,2*(k-1)+1],2),",",
                               round(con_result[[i]][3,2*(k-1)+1],2),")")
                out[i,2*k] <- 
                        paste0(round(con_result[[i]][1,2*k],1),
                               "(",round(con_result[[i]][2,2*k],1),",",
                               round(con_result[[i]][3,2*k],1),")")
        }
}

write.csv(out[,c(1,3,5)],paste0(path,"table_s6_update.csv"),row.names = F)

##
#####

## end of script

#####