#------------
# bootstrap Ct distribution
# GAM and skewness
# BY Lin Y. and Yang B.
# August 2021
#------------
######################################################
## data_ct: all individual Ct values (with test dates)
## data_daily_all: daily case counts/sample counts, incidence-based Rt; 
##                 daily Ct mean, median and skewness (imputed)
##                 correspond to "Supplementary" data in source data file
######################################################
# load packages
require(e1071) # to calculate skewness
require(mgcv)
#
# read in "data_ct.csv" and "data_daily_all.csv"
ct.linelist <- read.csv("data_ct.csv")
daily.linelist <- read.csv("data_daily_all.csv",as.is=T)
#
# function for bootstrap sampling
bootsample <- function (m,n,original.data){ 
        # randomly sampled M times, N samples each time
        mat.tmp <- matrix(NA,n,m) 
        #
        for (i in 1:m){
                newdata <- # randomly take N sample value from all records on that date
                        original.data[sample(1:nrow(original.data), n, 
                                             replace=TRUE), ]
                # store sampled values
                mat.tmp[,i] <- newdata$ct.value
        }
        #
        return(mat.tmp)
}
##
# set up time intervals 
start.date <- c("2020-07-01","2020-11-01")
end.date <- c("2020-08-31","2021-03-31")
#
## calculate GAM (for two waves separately)
# for daily GAM
wave.list <- list() 
# for bootstrap
ct.bs.wave <- ct.bs.daily <- list()
set.seed(615) # set seed to get same results each time
for (i in 1:2){
        date.seq <- seq(as.Date(start.date[i]),as.Date(end.date[i]),1)
        ct.tmp <- # all Ct records during each study waves
                ct.linelist[ct.linelist$date.test%in%as.character(date.seq),]
        ct.tmp$date.num <- as.numeric(as.Date(ct.tmp$date.test)) # time indicator for smoothing
        gam.used <- gam(ct.value~s(date.num),data=ct.tmp)
        daily.tmp <- daily.linelist[daily.linelist$date%in%as.character(date.seq),]
        daily.tmp$ct.gam <- predict(gam.used,daily.tmp)
        wave.list[[i]] <- daily.tmp 
        # take daily bootstrap samples 
        #(daily sample size equal to actual daily number of sample collection)
        for (ii in 1:length(date.seq)){
                ct.daily <- 
                        ct.tmp[ct.tmp$date.test==as.character(date.seq[ii]),]
                if (nrow(ct.daily)!=0){
                        n.take <- nrow(ct.daily)
                        ct.bs.daily[[ii]] <- bootsample(m=500,n=n.take,
                                                        original.data=ct.daily) 
                } else {
                        ct.bs.daily[[ii]] <- NA
                }
        }
        # store bootstrap records by waves
        ct.bs.wave[[i]] <- ct.bs.daily
}
#
## store bootstrap results in the data frame
replicate.n <- max(daily.linelist$records) + 10 
#
df.bs.list <- list()
for (i in 1:2){
        ct.bs.list <- ct.bs.wave[[i]]
        date.seq <- seq(as.Date(start.date[i]),as.Date(end.date[i]),1)
        df.bootstrap <- data.frame(date=rep(date.seq,each=replicate.n))
        df.bootstrap$date.num <- as.numeric(as.Date(df.bootstrap$date))
        df.bootstrap[,paste0("ct.value.",1:500)] <- NA
        for (n in 1:length(date.seq)){ # n days
                ct.num.tmp <- ct.bs.list[[n]]
                for (ii in 1:500){ # bootstrap 500 times
                        if (!all(is.na(ct.num.tmp))){
                                # one column per time
                                df.bootstrap[,paste0("ct.value.",ii)][(replicate.n*(n-1)+1):
                                                                             (replicate.n*(n-1)+nrow(ct.num.tmp))] <- ct.num.tmp[,ii] 
                                df.bootstrap[,paste0("ct.value.",ii)][(replicate.n*(n-1)+nrow(ct.num.tmp)+1):
                                                                             (replicate.n*n)] <- NA
                        } else {
                                df.bootstrap[,paste0("ct.value.",ii)][(replicate.n*(n-1)+1):(replicate.n*n)] <- NA
                        }
                }
        }
        df.bootstrap <- df.bootstrap[!is.na(df.bootstrap$ct.value.1),] # remove NAs
        df.bs.list[[i]] <- df.bootstrap
}
#
## get CI for GAM Ct and skewness
value.bs.list <- list()
for (i in 1:2){
        df.bootstrap <- df.bs.list[[i]]
        all.ct.bs <- data.frame(date=seq(as.Date(start.date[i]),
                                         as.Date(end.date[i]),1))
        all.ct.bs$date.num <- as.numeric(as.Date(all.ct.bs$date))
        all.ct.bs[,paste0("skewness.",1:500)] <- 
                all.ct.bs[,paste0("gam.",1:500)] <- NA
        #
        for (ii in 1:500){
                ct.used <- df.bootstrap
                colnames(ct.used)[ii+2] <- "ct.value" # need to rename
                # change the column name of the used Ct value for GAM
                gam.tmp <- gam(ct.value~s(date.num), data = ct.used)
                # get the GAM Ct 
                all.ct.bs[,paste0("gam.",ii)] <- predict(gam.tmp,all.ct.bs)
                ##
                ###
                for (n in 1:nrow(all.ct.bs)){
                        ct.tmp <- ct.used[ct.used$date==as.character(all.ct.bs$date[n]),]
                        if (nrow(ct.tmp)!=0){
                                all.ct.bs[,paste0("skewness.",ii)][n] <- 
                                        e1071::skewness(ct.tmp$ct.value)
                        }
                }
        }
        #
        value.bs.list[[i]] <- all.ct.bs
}
#
## get corresponding percentile CIs and merge with wave.list
wave.list.update <- list()
for (i in 1:2){
        all.out <- wave.list[[i]]
        all.bs.tmp <- value.bs.list[[i]]
        all.out$skewness.ub <- all.out$skewness.lb <- 
                all.out$gam.ub <- all.out$gam.lb <- NA # output
        for (n in 1:nrow(all.out)){
                ct.grab <- 
                        as.numeric(all.bs.tmp[,paste0("gam.",1:500)][
                                all.bs.tmp$date==as.character(all.out$date[n]),])
                # 95% confidence interval
                all.out$gam.lb[n] <- quantile(ct.grab,.025)
                all.out$gam.ub[n] <- quantile(ct.grab,.975)
                #        
                skewness.grab <- 
                        as.numeric(all.bs.tmp[,paste0("skewness.",1:500)][
                                all.bs.tmp$date==as.character(all.out$date[n]),])
                # 95% confidence interval
                all.out$skewness.lb[n] <- quantile(skewness.grab,.025,na.rm = T)
                all.out$skewness.ub[n] <- quantile(skewness.grab,.975,na.rm = T)
        }
        all.out$test.to.start <- as.numeric(as.Date(all.out$date)-
                                                    as.Date(start.date[i]))
        wave.list.update[[i]] <- all.out
}
#
wave.all <- rbind(wave.list.update[[1]],wave.list.update[[2]])
wave.all$period <- ifelse(as.Date(wave.all$date)>=start.date[1]&
                                  as.Date(wave.all$date)<=end.date[1],1,2)
#
## export results
# wave.all - with CIs calculated for GAM Ct and skewness
#write.csv(wave.all,"daily_ct_bootstrap.csv",row.names=F)
# corresponded to "Figure 1" data in source data file
##
#####

## end of script

#####
