#------------
# log-linear
# Rt and Ct
# also include reverse validation
#------------
######################################################
## data_daily_all: daily case counts/sample counts, incidence-based Rt; 
##                 daily Ct mean, median and skewness (imputed) from "1_merge_data"
######################################################
# read in "data_daily_all.csv"
daily.linelist <- read.csv("data_daily_all.csv",as.is=T)
#
######
# correlation
#####
# assign training/testing periods --
daily.linelist$period <- 
        ifelse(as.Date(daily.linelist$date)>=as.Date("2020-07-01")&
                       as.Date(daily.linelist$date)<=as.Date("2020-08-31"),1,
               ifelse(as.Date(daily.linelist$date)>=as.Date("2020-11-01")&
                              as.Date(daily.linelist$date)<=as.Date("2021-03-31"),2,0))
table(daily.linelist$period,useNA = 'always') # checked
#
# function for calculating spearman rho --
correlation.rho <- function(df,var1,var2){
        cortest <- cor.test(df[,var1],df[,var2],
                            use="na.or.complete",method="spearman")
        out <- round(c(cortest$estimate,cortest$p.value),2)
        return(out)
}
# calculate rho between Ct mean/skewness and Rt
cor.mat <- matrix(NA,2,4)
for (i in 1:2){
        df.tmp <- daily.linelist[daily.linelist$period==i,]
        df.tmp$log.rt <- log(df.tmp$local.rt.mean)
        cor.mat[1,(2*i-1):(2*i)] <- 
                correlation.rho(df.tmp, var1 = "mean", var2 = "log.rt")
        cor.mat[2,(2*i-1):(2*i)] <- 
                correlation.rho(df.tmp, var1 = "skewness", var2 = "log.rt")
}
cor.mat[cor.mat==0] <- "<0.001"
## export results
# cor.mat - Supplementary Table 1
#write.csv(cor.mat,"table_s1.csv",row.names=F)
#--------------
#####
# regression
#####
train.period <- seq(as.Date("2020-07-06"),as.Date("2020-08-31"),1)
train.data <- daily.linelist[daily.linelist$date%in%as.character(train.period),]
#
# model select on AIC
# original form, Rt
aic.list <- list()
m1 <- lm(local.rt.mean~mean,data=train.data)
m2 <- lm(local.rt.mean~median,data=train.data)
m3 <- lm(local.rt.mean~skewness.imputed,data=train.data)
m4 <- lm(local.rt.mean~mean+skewness.imputed,data=train.data)
m5 <- lm(local.rt.mean~median+skewness.imputed,data=train.data)
aic.list[[1]] <- AIC(m1,m2,m3,m4,m5)
#
# log-scaled Rt
n1 <- lm(log(local.rt.mean)~mean,data=train.data)
n2 <- lm(log(local.rt.mean)~median,data=train.data)
n3 <- lm(log(local.rt.mean)~skewness.imputed,data=train.data)
n4 <- lm(log(local.rt.mean)~mean+skewness.imputed,data=train.data)
n5 <- lm(log(local.rt.mean)~median+skewness.imputed,data=train.data)
aic.list[[2]] <- AIC(n1,n2,n3,n4,n5)
#
aic.mat <- matrix(NA,5,2)
for (i in 1:5){
        for (k in 1:2){
                aic.mat[i,k] <- aic.list[[k]]$AIC[i]
        }
}
aic.mat <- round(aic.mat,2)
## export results
# aic.mat - Supplementary Table 4
#write.csv(aic.mat,"table_s4.csv",row.names=F)
#--------------
#####
# get estimate
#####
model.used <- lm(log(local.rt.mean)~mean+skewness.imputed,data=train.data)
round(exp(cbind(coef(model.used),confint(model.used))),2)
summary(model.used) # get model coefficients and adjusted R^2
#
# get daily Ct-based Rt (as main estimates)
est <- exp(predict(model.used,daily.linelist,interval = "prediction"))
daily.est <- cbind(daily.linelist,est)
#
# Spearman rank correlation coefficients (rho)
# between incidence- and Ct-based Rt
# 
cor.rt <- matrix(NA,2,2)
for (i in 1:2){
        cortest.tmp <- 
                with(daily.est[daily.est$period==i,],
                     cor.test(log(local.rt.mean),log(fit),
                              use="na.or.complete",method="spearman"))
        cor.rt[i,1] <- round(cortest.tmp$est,2)
        cor.rt[i,2] <- round(cortest.tmp$p.value,2)
}
cor.rt
#
## export estimated daily Ct-based Rt
# daily.est - daily estimated Ct-based Rt
#write.csv(daily.est,"daily_ct_rt.csv",row.names = F)
##
#--------------
#####
# get coefficients
# for regression
#####
#
# get models build over the alternative training period 
train.period2 <- seq(as.Date("2020-11-20"),as.Date("2020-12-31"),1)
train.data2 <- daily.linelist[daily.linelist$date%in%as.character(train.period2),]
model.validate <- lm(log(local.rt.mean)~mean+skewness.imputed,data=train.data2)
#
model.list <- list(model.used,model.validate)
coef.mat <- matrix(NA,4,4)
for (i in 1:2){
        # coefficients
        coef.mat[(2*i-1):(2*i),1:3] <- 
                round(cbind(exp(coef(model.list[[i]])[2:3]),
                            exp(confint(model.list[[i]])[2:3,1:2])),2)
        # p-values
        coef.mat[(2*i-1):(2*i),4] <- 
                round(summary(model.list[[i]])$coefficients[2:3,4],2)
}
coef.out <- matrix(NA,4,3)
for (i in 1:4){
        coef.out[i,1] <- paste0(coef.mat[i,1],"(",coef.mat[i,2],
                                ",",coef.mat[i,3],")")
        coef.out[i,2] <- coef.mat[i,4]
}
# add r-square
coef.out[1,3] <- round(summary(model.list[[1]])$adj.r.square,2)
coef.out[3,3] <- round(summary(model.list[[2]])$adj.r.square,2)
coef.out[is.na(coef.out)] <- ""
coef.out[coef.out==0] <- "<0.001"
## export results
# coef.out - Supplementary Table 2
#write.csv(coef.out,"table_s2.csv",row.names = F)
##
#####

## end of script

#####
