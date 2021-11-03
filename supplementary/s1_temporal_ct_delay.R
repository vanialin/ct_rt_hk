#------------
# temporal distribution of
# first Ct and testing delays
# Fig S1
# BY Lin Y.
# updated October 2021
#------------
#
# load packages
require(plotrix)
require(lubridate)
require(dplyr)
require(scales)
#
######################################################
## data_ct: all individual Ct values (with test dates)
## daily_ct_bootstrap: daily case counts/sample counts, incidence-based Rt; 
##                     daily Ct mean, median and skewness (imputed)
##                     CIs calculated for GAM Ct and skewness (from "2_ct_for_bootstrap")
######################################################
# read in "data_ct.csv" and "data_cases.csv"
ct.linelist <- read.csv("/Users/vanialam/OneDrive - connect.hku.hk/vanialam/research_vania/epi_wave_2021/program/2021_08_R0/publish/data/data_ct.csv")
daily.ct <- read.csv("/Users/vanialam/OneDrive - connect.hku.hk/vanialam/research_vania/epi_wave_2021/program/2021_08_R0/publish/result/daily_ct_bootstrap.csv")
#

# plotting elements
# elements
# x-axis
start.date <- c("2020-07-01","2020-11-01")
end.date <- c("2020-08-31","2021-03-31")
x.day.list <- x.month.pos <- x.month.list <- x.pos.list <- list()
x.length <- rep(NA,2)
month.end <- list(c("2020-07-31","2020-08-31"),
                  c("2020-11-30","2020-12-31","2021-01-31","2021-02-28","2021-03-31"))
for (i in 1:2){
        # start of each 7-day interval
        if (i == 1){
                x.tmp <- c(as.Date(start.date[i])+0:8*7)
        } else {
                x.tmp <- c(as.Date(start.date[i])+0:21*7)
        }
        x.seq <- seq(as.Date(start.date[i]),as.Date(end.date[i]),1)
        x.length[i] <- length(x.seq)
        x.pos.list[[i]] <- # location for ticks
                which(x.seq%in%x.tmp)-1 # 
        x.day.list[[i]] <- day(x.tmp) # text for labelling
        x.month.list[[i]] <- which(day(x.seq)==15)-2 # location for labelling
        x.month.pos[[i]] <- c(0,(which(as.character(x.seq)%in%month.end[[i]])-1)) # for axis tick
}
## adding the gap for wave 4
gap.wave4 <- x.length[1]+9
x.pos.list[[2]] <- x.pos.list[[2]]+gap.wave4
x.month.pos[[2]] <- x.month.pos[[2]]+gap.wave4
x.month.list[[2]] <- x.month.list[[2]]+gap.wave4
# combine
x.total <- x.length[1]+x.length[2]+9
x.pos.combine <- c(x.pos.list[[1]],x.pos.list[[2]])
x.day.combine <- c(x.day.list[[1]],x.day.list[[2]])
x.month.combine <- c(x.month.list[[1]],x.month.list[[2]])
x.month.pos.combine <- c(x.month.pos[[1]],x.month.pos[[2]])
#
x.month.lab <- list(c("July","Aug"),c("Nov","Dec","Jan","Feb","Mar"))
month.title <- c("Training period","Testing period")
start.vec <- c(0,71)
#

### prepare the median for both delay and Ct
median.ct.list <- median.delay.list <- mbs.list <- list()
for (i in 1:2){
        df.used <- # sort by date.test
                ct.linelist[as.numeric(as.Date(ct.linelist$date.test)-
                                                 as.Date(start.date[i]))>=0 & 
                                      as.numeric(as.Date(ct.linelist$date.test)-
                                                         as.Date(end.date[i]))<=0,]
        df.used$test.to.start <- 
                as.numeric(as.Date(df.used$date.test)-as.Date(start.date[i]))
        df.used <- df.used[df.used$test.to.onset<=40,] # for plotting
        mbs.list[[i]] <- df.used
        a <- boxplot(df.used$test.to.onset~df.used$test.to.start)
        median.delay.list[[i]] <- a$stats[3,]
        b <- boxplot(df.used$ct.value~df.used$test.to.start)
        median.ct.list[[i]] <- b$stats[3,]
}
#
##
## get GAM delay and CIs
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
                mat.tmp[,i] <- newdata$test.to.onset
        }
        #
        return(mat.tmp)
}
##
#
## calculate GAM (for two waves separately)
# for daily GAM
wave.list <- list() 
# for bootstrap
delay.bs.wave <- delay.bs.daily <- list()
set.seed(615) # set seed to get same results each time
for (i in 1:2){
        date.seq <- seq(as.Date(start.date[i]),as.Date(end.date[i]),1)
        delay.tmp <- # all Ct records during each study waves
                ct.linelist[ct.linelist$date.test%in%as.character(date.seq)&
                                    !is.na(ct.linelist$date.onset),]
        delay.tmp$date.num <- as.numeric(as.Date(delay.tmp$date.test)) # time indicator for smoothing
        gam.used <- gam(test.to.onset~s(date.num),data=delay.tmp)
        daily.tmp <- daily.ct[daily.ct$date%in%as.character(date.seq),]
        daily.tmp$delay.gam <- predict(gam.used,daily.tmp)
        wave.list[[i]] <- daily.tmp 
        # take daily bootstrap samples 
        #(daily sample size equal to actual daily number of sample collection)
        for (ii in 1:length(date.seq)){
                delay.daily <- 
                        delay.tmp[delay.tmp$date.test==as.character(date.seq[ii]),]
                if (nrow(delay.daily)!=0){
                        n.take <- nrow(delay.daily)
                        delay.bs.daily[[ii]] <- bootsample(m=500,n=n.take,
                                                        original.data=delay.daily) 
                } else {
                        delay.bs.daily[[ii]] <- NA
                }
        }
        # store bootstrap records by waves
        delay.bs.wave[[i]] <- delay.bs.daily
}
#
## store bootstrap results in the data frame
replicate.n <- max(daily.ct$records) + 10 
#
df.bs.list <- list()
for (i in 1:2){
        delay.bs.list <- delay.bs.wave[[i]]
        date.seq <- seq(as.Date(start.date[i]),as.Date(end.date[i]),1)
        df.bootstrap <- data.frame(date=rep(date.seq,each=replicate.n))
        df.bootstrap$date.num <- as.numeric(as.Date(df.bootstrap$date))
        df.bootstrap[,paste0("delay.",1:500)] <- NA
        for (n in 1:length(date.seq)){ # n days
                delay.num.tmp <- delay.bs.list[[n]]
                for (ii in 1:500){ # bootstrap 500 times
                        if (!all(is.na(delay.num.tmp))){
                                # one column per time
                                df.bootstrap[,paste0("delay.",ii)][(replicate.n*(n-1)+1):
                                                                              (replicate.n*(n-1)+nrow(delay.num.tmp))] <- delay.num.tmp[,ii] 
                                df.bootstrap[,paste0("delay.",ii)][(replicate.n*(n-1)+nrow(delay.num.tmp)+1):
                                                                              (replicate.n*n)] <- NA
                        } else {
                                df.bootstrap[,paste0("delay.",ii)][(replicate.n*(n-1)+1):(replicate.n*n)] <- NA
                        }
                }
        }
        df.bootstrap <- df.bootstrap[!is.na(df.bootstrap$delay.1),] # remove NAs
        df.bs.list[[i]] <- df.bootstrap
}
#
## get CI for GAM delay
value.bs.list <- list()
for (i in 1:2){
        df.bootstrap <- df.bs.list[[i]]
        all.delay.bs <- data.frame(date=seq(as.Date(start.date[i]),
                                         as.Date(end.date[i]),1))
        all.delay.bs$date.num <- as.numeric(as.Date(all.delay.bs$date))
        all.delay.bs[,paste0("gam.",1:500)] <- NA
        #
        for (ii in 1:500){
                delay.used <- df.bootstrap
                colnames(delay.used)[ii+2] <- "delay" # need to rename
                # change the column name of the used Ct value for GAM
                gam.tmp <- gam(delay~s(date.num), data = delay.used)
                # get the GAM Ct 
                all.delay.bs[,paste0("gam.",ii)] <- predict(gam.tmp,all.delay.bs)
                ##
        }
        #
        value.bs.list[[i]] <- all.delay.bs
}
#
## get corresponding percentile CIs and merge with wave.list
wave.list.update <- list()
for (i in 1:2){
        all.out <- wave.list[[i]]
        all.bs.tmp <- value.bs.list[[i]]
        all.out$delay.ub <- all.out$delay.lb <- NA # output
        for (n in 1:nrow(all.out)){
                delay.grab <- 
                        as.numeric(all.bs.tmp[,paste0("gam.",1:500)][
                                all.bs.tmp$date==as.character(all.out$date[n]),])
                # 95% confidence interval
                all.out$delay.lb[n] <- quantile(delay.grab,.025)
                all.out$delay.ub[n] <- quantile(delay.grab,.975)
        }
        all.out$test.to.start <- as.numeric(as.Date(all.out$date)-
                                                    as.Date(start.date[i]))
        wave.list.update[[i]] <- all.out
}
#

#
#

### plot (8*12)
#par(mar=c(4,3,2,0)+0.1)
par(fig=c(0,1,0.45,1))
plot(NA,xlim=c(1,x.total),ylim=rev(c(10,40)),xlab=NA,ylab=NA,axes = F,main=NA)
axis(2,at=2:8*5,las=1,line=-.5)
mtext("Ct value",side=2,line=1.5)
# x-axis
day.axis <- c(0:(x.length[1]-1),gap.wave4:(x.length[2]+gap.wave4-1))
axis(1,at=day.axis,labels = rep(NA,length(day.axis)),tck=-.007)
axis(1,at=c(0,131,221),labels=rep(NA,3),tck=-.05)
axis.break(1,x.length[1]+4,style="slash",brw=.015)
# add axis legend separately
for (i in 1:2){
        axis(1,at=x.month.pos[[i]],
             labels = rep(NA,length(x.month.pos[[i]])),tck=-.025)
}
# Ct boxplots 
for (i in 1:2){
        mbs.tmp <- mbs.list[[i]]
        mbs.tmp$test.to.start <- mbs.tmp$test.to.start+start.vec[i]
        ord <- unique(mbs.tmp$test.to.start)[order(unique(mbs.tmp$test.to.start))]%>%na.omit()
        boxplot(mbs.tmp$ct.value~mbs.tmp$test.to.start,
                add=T,axes=F,at=ord,cex.lab=1.3,border=NA,
                boxwex=.2,boxcol="dark grey",whisklty = 1,staplecol="white",
                whiskcol=alpha("grey",.45),whisklwd=1.35)
        #median
        median.stat <- median.ct.list[[i]]
        for (ii in 1:length(median.stat)){
                points(ord[ii],
                       median.stat[ii],col="#ff4d00",pch=19,cex=.35)
        }
        # GAM Ct
        df.tmp <- wave.list.update[[i]]
        polygon(c(df.tmp$test.to.start+start.vec[i],
                  rev(df.tmp$test.to.start)+start.vec[i]),
                c(df.tmp$gam.lb,rev(df.tmp$gam.ub)),
                col=alpha("orange",.3),border=F) 
        lines(df.tmp$test.to.start+start.vec[i],df.tmp$ct.gam,col="orange",lwd=1.5)
}

##

## delay panel

par(fig=c(0,1,0,0.55),new=T)
plot(NA,xlim=c(1,x.total),ylim=c(0,15),xlab=NA,ylab=NA,axes = F,main=NA)
axis(2,at=0:3*5,las=1,line=-.5)
mtext("Delay",side=2,line=1.5)
# x-axis
day.axis <- c(0:(x.length[1]-1),gap.wave4:(x.length[2]+gap.wave4-1))
axis(1,at=day.axis,
     labels = rep(NA,length(day.axis)),tck=-.007)
axis(1,at=x.pos.combine,labels = x.day.combine,las=1,tck=-.012,padj=-1.2)
axis(1,at=c(0,131,221),labels=rep(NA,3),tck=-.05)
axis.break(1,x.length[1]+4,style="slash",brw=.015)
# add axis legend separately
for (i in 1:2){
        axis(1,at=x.month.pos[[i]],
             labels = rep(NA,length(x.month.pos[[i]])),tck=-.025)
        for (k in 1:length(x.month.list[[i]])){
                mtext(x.month.lab[[i]][k],side=1,line=1.2,
                      at=x.month.list[[i]][k],adj=0)
        } 
}
# add longer axis
mtext("2020",at = 66, side=1,line=2.5,font=2)
mtext("2021",at = 175, side=1,line=2.5,font=2)
# Ct boxplots 
for (i in 1:2){
        mbs.tmp <- mbs.list[[i]]
        mbs.tmp$test.to.start <- mbs.tmp$test.to.start+start.vec[i]
        ord <- unique(mbs.tmp$test.to.start)[order(unique(mbs.tmp$test.to.start))]%>%na.omit()
        boxplot(mbs.tmp$test.to.onset~mbs.tmp$test.to.start,
                add=T,axes=F,at=ord,cex.lab=1.3,border=NA,
                boxwex=.2,boxcol="dark grey",whisklty = 1,staplecol="white",
                whiskcol=alpha("grey",.45),whisklwd=1.35)
        #median
        median.stat <- median.delay.list[[i]]
        for (ii in 1:length(median.stat)){
                points(ord[ii],
                       median.stat[ii],col="dark green",pch=19,cex=.35)
        }
        # GAM delay
        df.tmp <- wave.list.update[[i]]
        polygon(c(df.tmp$test.to.start+start.vec[i],
                  rev(df.tmp$test.to.start)+start.vec[i]),
                c(df.tmp$delay.lb,rev(df.tmp$delay.ub)),
                col=alpha(1,.3),border=F) 
        lines(df.tmp$test.to.start+start.vec[i],df.tmp$delay.gam,lwd=1.5)
}

###
#####

## end of script

#####