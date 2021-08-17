#-----------
# merge daily data
# case count and Ct distribution
#-----------
######################################################
## data_ct: all individual Ct values (with test dates)
## data_cases: daily case counts/sample counts and incidence-based Rt
######################################################
# read in "data_ct.csv" and "data_cases.csv"
#setwd()
ct.linelist <- read.csv("data_ct.csv",as.is=T)
daily.linelist <- read.csv("data_cases.csv")
#
# calculate Ct parameters (mean, median, skewness)
ct <- ct.linelist %>% 
        group_by(date.test) %>%
        summarise(mean = mean(ct.value, na.rm = TRUE),
                  median = median(ct.value, na.rm = TRUE),
                  skewness = e1071::skewness(ct.value,
                                             na.rm = TRUE)) %>%
        ungroup()
#
ct$date.test <- as.Date(ct$date.test)
ct1 = ct %>%
        bind_cols(
                # impute
                skewness.imputed = sapply(ct$date.test, function(d){
                        
                        mean(ct$skewness[ct$date.test %in% c((d-7):(d-1))], 
                             na.rm = TRUE)
                        
                })
        )%>%
        mutate(date.num = as.numeric(date.test),
               # categorical Ct parameters
               mean.cat = cut(mean,
                              c(10, 20, 22, 24, 26, 40),
                              labels = 1:5),
               median.cat = cut(median,
                                c(10, 20, 22, 24, 26, 40),
                                labels = 1:5),
               skewness.cat = cut(skewness,
                                  c(-1, -0.3, 0, 0.3, 1.5),
                                  labels = 1:4)
        )

ct1$skewness.imputed[!is.na(ct1$skewness)] <- ct1$skewness[!is.na(ct1$skewness)]
#
# merge
daily.out <- merge(daily.linelist,ct1,by.x = "date",all.x = T,by.y = "date.test")
daily.out$date.num <- as.numeric(daily.out$date)
## export
# daily.out - all daily information required for regression 
# i.e., daily incidence-based Rt, daily Ct measured by mean/median and skewness (imputed)
#write.csv(daily.out,"DIRECTORY/data_daily_all.csv",row.names = F)
##
#####

## end of script

#####