#------------
# back-projecting Ct
# on illness onset
# Fig. S2-4
# By Lin.Y and Yang.B
# updated October 2021
#------------
######################################################
## data_ct: all individual Ct values (with test dates)
## data_cases: daily case counts/sample counts and incidence-based Rt
######################################################
#
# load packages
require(scales)
require(mgcv)
require(readxl)
require(dplyr)
#
#setwd("/Users/vanialam/OneDrive - connect.hku.hk/vanialam/research_vania/epi_wave_2021/program/2021_09_R1/publish (EDIT HERE)/")
# read in "data_ct.csv" and "data_cases.csv"
ct.linelist <- read.csv("data/data_ct.csv")
daily.linelist <- read.csv("data/data_cases.csv")
#

### extrapolate Ct value at onset ----
## log-linear regression
ct.symp <- ct.linelist[!is.na(ct.linelist$date.onset),]
lm1 <- lm(log(ct.value)~test.to.onset*age.gp+factor(key),data=ct.symp)
# log-linear trend of Ct over time
# key - individual indicator (i.e., patient ID)

summary(lm1)$adj.r.square # R2 = 0.60
df <- data.frame(est=coef(lm1))
s <- confint(lm1)
df$lwr <- s[,1]
df$upr <- s[,2]

coef.df <- data.frame(key.factor=row.names(df),value=as.numeric(df$est))
coef.df$ct.onset <- coef.df$age.gp <- coef.df$key <- NA
for (i in 4:nrow(coef.df)){
        coef.df$key[i] <- substr(as.character(coef.df$key.factor[i]),12,16)
        coef.df$age.gp[i] <- ct.symp$age.gp[ct.symp$key==coef.df$key[i]][1]
}
# to estimate the Ct at onset for each symptomatic case
coef.out <- coef.df[4:(nrow(coef.df)-1),]
tail(coef.out) # check
coef.out$test.to.onset <- 0
pred.value <- exp(predict(lm1,coef.out,interval="confidence"))
coef.out$ct.onset <- pred.value[,1]
coef.out$lwr <- pred.value[,2]
coef.out$upr <- pred.value[,3]
coef.out$ci <- coef.out$upr-coef.out$lwr

coef.out <- coef.out[,c("key","age.gp","ct.onset","lwr","upr","ci")]
#write.csv(coef.out,"results/ct_onset_singleobs.csv",row.names = F)

ct.linelist2 <- merge(ct.linelist,coef.out[,c("key","ct.onset","lwr","upr")],
                      by.x = "key",all.x = T,by.y = "key")

### provide temporal Ct trend for each age group (Fig. S2) ----
lm.gen <- lm(log(ct.value)~test.to.onset*age.gp,data=ct.symp)
cbind(exp(coef(lm.gen)),exp(confint(lm.gen)))
col.point <- c("dark grey","#5bc0de","#d9534f")
col.line <- c("black","#0057e7","red")
text.add <- c("0-18 yo","19-64 yo","65+ yo")
### plot
pdf("results/Fig_S2.pdf",height = 6, width = 8)
plot(NA,xlim=c(0,20),ylim=rev(c(7,40)),las=1,xlab="Days since onset",ylab="Ct value",axes=F)
axis(1,at=0:4*5)
axis(2,at=2:8*5,las=1)
for (i in 1:3){
        df.tmp <- ct.linelist2[ct.linelist2$age.gp==i,] 
        points(df.tmp$test.to.onset-0.2+0.2*(i-1),df.tmp$ct.value,
               pch=19,col=alpha(col.point[i],.2))
}
pred.df <- data.frame(test.to.onset=rep(0:20,3),age.gp=rep(1:3,each=21))
pred.df[,3:5] <- exp(predict(lm.gen,pred.df,interval = "confidence"))
for (i in 1:3){
        polygon(c(0:20,rev(0:20)),
                c(pred.df[1:21+21*(i-1),4],
                  rev(pred.df[1:21+21*(i-1),5])),col=alpha(col.line[i],.3),border=F)
        lines(0:20,pred.df[1:21+21*(i-1),3],col=col.line[i],lwd=2)
        ## legend
        lines(c(13,15),rep(7+1.5*i,2),col=col.line[i],lwd=2)
        points(12.5,7+1.5*i,col=alpha(col.point[i],.6),pch=19)
        text(15.5,7+1.5*i,text.add[i],adj=0)
}
dev.off()
#
#

### demonstration of differences in temporal trend of Ct at sampling and at onset ###
## for more clear visual comparison
ct.linelist2$date.onset.num <- as.numeric(as.Date(ct.linelist2$date.onset))
ct.linelist2$date.test.num <- as.numeric(as.Date(ct.linelist2$date.test))
count.training <- daily.linelist[month(daily.linelist$date)%in%7:8,]
# compare when records were sufficient
(d <- range(as.Date(as.character(count.training$date[count.training$records>=30])))) 
date.seq <- seq(as.Date(d[1]),as.Date(d[2]),1)
date.seq.num <- as.numeric(date.seq)

ct.tmp <- ct.linelist2[as.Date(ct.linelist2$date.test)%in%date.seq,]

count.tmp <- ct.linelist2 %>% group_by(date.test.num) %>% 
        summarise(all.records=n()) %>% ungroup() %>%
        rename(date.onset.num=date.test.num) %>% # to facilitate binding
        right_join(ct.linelist2 %>% group_by(date.onset.num) %>% 
                           summarise(symp.records=n()) %>% ungroup()) %>%
        rename(date.num=date.onset.num) %>% 
        # comparison only made for time period with sufficient record counts
        filter(date.num%in%date.seq.num) 

ct.df.actual <- ct.tmp
ct.df.onset <- ct.tmp[ct.tmp$date.onset.num%in%date.seq.num&
                                  !is.na(ct.tmp$date.onset.num),]

## compare regression coefficient against time 
#lm.onset <- lm(date.onset.num~ct.onset,data=ct.df.onset)
#lm.actual <- lm(date.test.num~ct.value,data=ct.df.actual)
#round(c(coef(lm.onset)[2],confint(lm.onset)[2,],summary(lm.onset)$coefficients[2,4]),3)
# p=0.43, insignificant trend for Ct at onset
#round(c(coef(lm.actual)[2],confint(lm.actual)[2,],summary(lm.actual)$coefficients[2,4]),3)
# p<0.001, significant trend for Ct at sampling

## prepare median for both distribution
a <- boxplot(ct.value~date.test.num,data=ct.df.actual)
median.actual <- a$stats[3,]
b <- boxplot(ct.onset~date.onset.num,data=ct.df.onset)
median.onset <- b$stats[3,]

###
### Fig. S3 ----
pdf("results/Fig_S3.pdf",height = 7, width = 12)
par(fig=c(0,0.67,0,1),mar=c(4,3,2,3)+0.1)
plot(NA,xlim=range(date.seq.num),ylim=c(0,150),xlab=NA,ylab=NA,axes=F)
axis(4,las=1,at=0:5*30)
mtext("Number of cases",side=4,line=2)
for (j in 1:nrow(count.tmp)){
        polygon(c(rep(count.tmp$date.num[j]-0.5,2),
                  rep(count.tmp$date.num[j]+0.5,2)),
                c(0,rep(count.tmp$all.records[j],2),0),
                col="#838584", # all records **BY SAMPLING**
                border="white")
        polygon(c(rep(count.tmp$date.num[j]-0.5,2),
                  rep(count.tmp$date.num[j]+0.5,2)),
                c(0,rep(count.tmp$symp.records[j],2),0),
                col=alpha("light blue",.75), # symp record **BY ONSET**
                border="white")
}
par(new=T)
plot(NA,xlim=range(date.seq.num),ylim=rev(c(18,28)),las=1,xlab=NA,
     ylab=NA,axes=F)
axis(2,las=1)
mtext("Ct value",side=2,line=2)
date.pos <- c(date.seq.num[1],date.seq.num[(1:floor(length(date.seq.num)/10))*10],
              date.seq[length(date.seq.num)])
date.label <- as.character(as.Date(date.pos,
                                   origin = "1970-01-01"))
axis(1,at=date.seq.num,labels=rep(NA,length(date.seq.num)),tck=-.015)
axis(1,at=date.pos,
     labels = paste0(day(date.label),"/",month(date.label)))
mtext("Date",side=1,line=2)
#
## actual
gam.actual <- gam(ct.value~s(date.test.num),data=ct.df.actual)
p1 <- predict(gam.actual, data.frame(date.test.num=date.seq.num), 
              type = "link", se.fit = TRUE)
est <- p1$fit
upr <- p1$fit + (2 * p1$se.fit) 
lwr <- p1$fit - (2 * p1$se.fit) 
lines(range(date.seq.num),rep(mean(ct.df.actual$ct.value),2),
      col=alpha("orange",.7),lwd=2,lty=2) # reference for comparison
polygon(c(date.seq.num,rev(date.seq.num)),
        c(lwr,rev(upr)),col=alpha("orange",.3),border=F)
lines(date.seq.num,est,col="orange",lwd=2)
## estimated at onset
gam.onset <- gam(ct.onset~s(date.onset.num),data=ct.df.onset)
p2 <- predict(gam.onset, data.frame(date.onset.num=date.seq.num), 
              type = "link", se.fit = TRUE)
est2 <- p2$fit
upr2 <- p2$fit + (2 * p2$se.fit) 
lwr2 <- p2$fit - (2 * p2$se.fit) 
lines(range(date.seq.num),rep(mean(ct.df.onset$ct.onset),2),
      col=alpha("blue",.7),lwd=2,lty=2) # reference for comparison
polygon(c(date.seq.num,rev(date.seq.num)),
        c(lwr2,rev(upr2)),col=alpha("blue",.3),border=F)
lines(date.seq.num,est2,col="blue",lwd=2)
## legend
polygon(c(rep(18477,2),rep(18478,2)),c(17.8,18,18,17.8),col="grey",border = F)
polygon(c(rep(18477,2),rep(18478,2)),c(18.1,18.3,18.3,18.1),col="light blue",border = F)
text(18478.5,17.9,"Case by sampling",adj=0)
text(18478.5,18.2,"Case by onset",adj=0)
lines(c(18477,18478),rep(18.5,2),col="orange",lwd=1.5)
lines(c(18477,18478),rep(18.8,2),col="blue",lwd=1.5)
text(18478.5,18.5,"GAM Ct by sampling",adj=0)
text(18478.5,18.8,"GAM Ct by onset",adj=0)
lines(c(18477,18478),rep(19.1,2),col="orange",lwd=1.5,lty=2)
lines(c(18477,18478),rep(19.4,2),col="blue",lwd=1.5,lty=2)
text(18478.5,19.1,"Mean value of Ct by sampling",adj=0)
text(18478.5,19.4,"Mean value of Ct by onset",adj=0)
mtext("a",side=3,font=2,adj=0,cex=1.5)
#
#
par(fig=c(0.67,1,0.45,1),new=T,mar=c(4,4,2,2)+0.1)
plot(NA,xlim=range(date.seq.num),ylim=rev(c(10,40)),xlab=NA,ylab=NA,axes = F,main=NA)
axis(2,las=1)
mtext("Ct value",side=2,line=2)
date.pos <- c(date.seq.num[1],date.seq.num[(1:floor(length(date.seq.num)/10))*10],
              date.seq.num[length(date.seq.num)])
date.label <- as.character(as.Date(date.pos,
                                   origin = "1970-01-01"))
axis(1,at=date.seq.num,labels=rep(NA,length(date.seq.num)),tck=-.015)
axis(1,at=date.pos,labels = rep(NA,length(date.label)))

ord <- unique(ct.df.actual$date.test.num)[order(unique(ct.df.actual$date.test.num))]
boxplot(ct.df.actual$ct.value~ct.df.actual$date.test.num,
        add=T,axes=F,at=ord,cex.lab=1.3,border=NA,
        boxwex=.2,boxcol="dark grey",whisklty = 1,staplecol="white",
        whiskcol=alpha("grey",.45),whisklwd=1.65)
for (ii in 1:length(median.actual)){
        points(ord[ii],median.actual[ii],col="orange",pch=19,cex=.45)
}
polygon(c(date.seq.num,rev(date.seq.num)),
        c(lwr,rev(upr)),col=alpha("orange",.3),border=F)
lines(date.seq.num,est,col="orange",lwd=2)
mtext("b",side=3,font=2,adj=0,cex=1.5)
#
#
par(fig=c(0.67,1,0,0.55),new=T)
plot(NA,xlim=range(date.seq.num),ylim=rev(c(10,40)),xlab=NA,ylab=NA,axes = F,main=NA)
axis(2,las=1)
mtext("Ct value",side=2,line=2)
date.pos <- c(date.seq.num[1],date.seq.num[(1:floor(length(date.seq.num)/10))*10],
              date.seq.num[length(date.seq.num)])
date.label <- as.character(as.Date(date.pos,
                                   origin = "1970-01-01"))
axis(1,at=date.seq.num,labels=rep(NA,length(date.seq.num)),tck=-.015)
axis(1,at=date.pos,
     labels = paste0(day(date.label),"/",month(date.label)))
mtext("Date",side=1,line=2)
ord2 <- unique(ct.df.onset$date.onset.num)[order(unique(ct.df.onset$date.onset.num))]
boxplot(ct.df.onset$ct.onset~ct.df.onset$date.onset.num,
        add=T,axes=F,at=ord2,cex.lab=1.3,border=NA,
        boxwex=.2,boxcol="light blue",whisklty = 1,staplecol="white",
        whiskcol=alpha("light blue",.45),whisklwd=1.65)
for (ii in 1:length(median.onset)){
        points(ord2[ii],median.onset[ii],col="blue",pch=19,cex=.35)
}
polygon(c(date.seq.num,rev(date.seq.num)),
        c(lwr2,rev(upr2)),col=alpha("blue",.3),border=F)
lines(date.seq.num,est2,col="blue",lwd=1.5)
mtext("c",side=3,font=2,adj=0,cex=1.5)
dev.off()

#
#

### the whole study period ----
## compare skewness
#
start.date <- c("2020-07-01","2020-11-01")
end.date <- c("2020-08-31","2021-03-31")
#
ct.wave <- list() # store by waves
for (i in 1:2){ # sort by date.test
        date.seq <- seq(as.Date(start.date[i]),as.Date(end.date[i]),1)
        ct.wave[[i]] <- ct.linelist2[ct.linelist2$date.test%in%
                                             as.character(date.seq),]
}
# prepare plotting elements
cut.date1 <- seq(as.Date(start.date[1]),as.Date(end.date[1]),14)
cut.date2 <- seq(as.Date(start.date[2]),as.Date(end.date[2]),14)
cut.date.list <- list(cut.date1,cut.date2)
actual.list <- onset.list <- list()
avg.rt <- skew.actual <- skew.onset <- rep(NA,14)
title.start <- title.end <- rep(NA,14) # for writing the title for histogram
for (i in 1:2){
        ct.tmp <- ct.wave[[i]]
        cut.date <- cut.date.list[[i]]
        for (n in 1:(length(cut.date)-1)){
                date.seq <- seq(as.Date(cut.date[n]),as.Date(cut.date[n+1]),1)
                actual.used <- 
                        ct.tmp[as.character(ct.tmp$date.test)%in%
                                       as.character(date.seq),]
                onset.used <-  # by **DATE.ONSET** 
                        ct.tmp[as.character(ct.tmp$date.onset)%in%
                                       as.character(date.seq)&
                                       !is.na(ct.tmp$date.onset),]
                # to store
                actual.list[[4*(i-1)+n]] <- actual.used
                onset.list[[4*(i-1)+n]] <- onset.used
                # calculate the skewness
                skew.actual[4*(i-1)+n] <- e1071::skewness(actual.used$ct.value)
                skew.onset[4*(i-1)+n] <- e1071::skewness(onset.used$ct.onset,na.rm = T)
                ## store date
                title.start[4*(i-1)+n] <- as.character(cut.date[n])
                title.end[4*(i-1)+n] <- as.character(as.Date(cut.date[n+1])-1)
                #
                #
                # to check epidemic situation
                data.tmp <- daily.linelist[as.Date(daily.linelist$date)%in%date.seq,]
                diff.vec <- diff(data.tmp$local.rt.mean)
                if (sum(diff.vec<0)>=(nrow(data.tmp)/2)&
                    min(data.tmp$local.rt.mean)<=1){ # declining
                        avg.rt[4*(i-1)+n] <- 1 # low epidemic
                } else {
                        avg.rt[4*(i-1)+n] <- 0 # high epidemic
                }
        }
}
skew.actual;skew.onset
summary(skew.actual);IQR(skew.actual)
summary(skew.onset);IQR(skew.onset)
#
#
### Fig. S4 ----
pdf("results/Fig_S4.pdf",height = 10,width = 9)
par(mar=c(4,3,2,2)+0.1)
x.vec <- c(rep(0.4,4),rep(0.7,3),rep(1,3))
y.vec <- c(c(1,0.77,0.54,0.31),rep(c(1,0.77,0.54),2))
col.option <-  c("#eea990","orange")
col.option2 <- c("#008080","#008744")
take <- c(1:4,6:8,9:11)
for (i in 1:10){
        if (i == 1){
                par(fig=c(x.vec[i]-0.37,x.vec[i],y.vec[i]-0.28,y.vec[i]))
        } else {
                par(fig=c(x.vec[i]-0.37,x.vec[i],y.vec[i]-0.28,y.vec[i]),new=T)
        }
        take.tmp <- take[i]
        # set threshold for distinguishing color (for actual Ct)
        # choose color
        if (avg.rt[take.tmp]==1){ # low epidemic
                col.actual <- col.option[1]
                col.onset <- col.option2[1]
        } else {
                col.actual <- col.option[2]
                col.onset <- col.option2[2]
        }
        #
        hist(actual.list[[take.tmp]]$ct.value,breaks = 4:20*2,
             col=alpha(col.actual,.95),border="white",main=NA,xlab=NA,ylab=NA,las=1,
             ylim = c(0,200))
        hist(onset.list[[take.tmp]]$ct.onset[onset.list[[take.tmp]]$ct.onset<=40],
             breaks = 4:20*2,col=alpha(col.onset,.55+0.15*(avg.rt[take.tmp]==1)),
             border="white",add=T,main=NA,ylim = c(0,200))
        abline(v=mean(actual.list[[take.tmp]]$ct.value),col="#d46c26",lwd=2)
        abline(v=mean(onset.list[[take.tmp]]$ct.onset,na.rm = T),
               col="#008744",lwd=2,lty=2)
        # legends (for skewness)
        # actual
        mtext(expression(b[s]),side=3,line=-1.2,at=30,col=col.actual,cex=.8)
        mtext(paste0("= ",round(skew.actual[take.tmp],2)),side=3,
              line=-1.2,at=31.5,col=col.actual,adj=0,cex=.8)
        # onset
        mtext(expression(b[o]),side=3,line=-2.2,at=30,col=col.onset,cex=.8)
        mtext(paste0("= ",round(skew.onset[take.tmp],2)),side=3,
              line=-2.2,at=31.5,col=col.onset,adj=0,cex=.8)
        mtext(paste0(day(title.start[take.tmp]),"/",
                     month(title.start[take.tmp]),"-",
                     day(title.end[take.tmp]),"/",
                     month(title.end[take.tmp])),side=3,line=0,cex=.9,font=2)
        # axis legend
        if (i %in% 1:4){
                mtext("Number",side=2,line=2.4,cex=.9)
        }
        if (i %in% c(4,7,10)){
                mtext("Ct",side=1,line=2,cex=.9)
        }
        # legend
        if (i == 1){
                mtext("a",side=3,adj=0,font=2,line=.5,at=.9,cex=1.3)
        } else if (i == 5){
                mtext("b",side=3,adj=0,font=2,line=.5,at=.9,cex=1.3)
        } else if (i == 8){
                mtext("c",side=3,adj=0,font=2,line=.5,at=.9,cex=1.3)
        }
        
}
# legend
par(fig=c(0.45,1,0,0.31),new=T)
plot(NA,xlim=c(0,10),ylim=c(0,10),axes=F,xlab=NA,ylab=NA)
# increasing
text(0,9.8,"Increasing epidemic",adj=0)
polygon(c(rep(0,2),rep(1,2)),c(8.8,9,9,8.8),col="orange", # all record by sampling
        border="white")
text(1.3,8.9,"Ct at sampling",adj=0)
polygon(c(rep(0,2),rep(1,2)),c(7.8,8,8,7.8),col=alpha(col.option2[2],.55), # all record by onset
        border="white")
text(1.3,7.9,"Ct at onset",adj=0)
###
# decreasing
text(0,6.5,"Decreasing epidemic",adj=0)
polygon(c(rep(0,2),rep(1,2)),c(5.5,5.7,5.7,5.5),col=col.option[1], # all record by sampling
        border="white")
text(1.3,5.6,"Ct at sampling",adj=0)
polygon(c(rep(0,2),rep(1,2)),c(4.5,4.7,4.7,4.5),col=alpha(col.option2[1],.7), # all record by onset
        border="white")
text(1.3,4.6,"Ct at onset",adj=0)
dev.off()
#
#
# calculate coefficients of variation
sd(skew.actual)/mean(skew.actual)
sd(skew.onset)/mean(skew.onset)
##
#####

## end of script

#####