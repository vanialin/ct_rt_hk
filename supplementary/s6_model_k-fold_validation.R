#------------
# model validation
# 10-fold cross validation
# By Yang.B and Lin.Y
# August 2021
#------------
######################################################
## data_daily_all: daily case counts/sample counts, incidence-based Rt; 
##                 daily Ct mean, median and skewness (imputed)
######################################################
# read in "data_daily_all.csv"
#daily.linelist <- read.csv("/Users/vanialam/OneDrive - connect.hku.hk/vanialam/research_vania/epi_wave_2021/program/publish/data/data_daily_all.csv",as.is=T)
#
# data on days for cross-validation
data1 <- daily.linelist[as.Date(daily.linelist$date)>=
                                as.Date("2020-07-06")&
                                # only validate on days with no less than 5 samples
                                daily.linelist$records>5,]
#
# 10-fold validation rounds
set.seed(730) # set seed to get same results every time
split1 <- split(data1, sample(1:nrow(data1), 10, replace=F)) # split into 10 validation sets
validate.mae <- consistency <- rep(NA,10)
for (i in 1:10){
        validate.set <- split1[[i]]
        # trained on days not covered in that validation round
        train.set <- data1[!data1$date%in%(split1[[i]]$date),]
        # main model --
        lm.tmp <- lm(log(local.rt.mean)~mean+skewness.imputed,data=train.set)
        validate.set$fit <- exp(predict(lm.tmp,validate.set))
        # calculate MAE
        validate.mae[i] <- 
                sum(abs(log(validate.set$fit)-
                                log(validate.set$local.rt.mean)),na.rm = T)/
                sum(!is.na(validate.set$fit))
        # calculate consistency
        consistency[i] <- 
                nrow(validate.set[((validate.set$local.rt.mean-1)*(validate.set$fit-1))>0,])/
                nrow(validate.set)
}
#
mean(consistency);range(consistency)
validate.mae
mean(validate.mae);range(validate.mae) # 
##
#####

## end of script

#####