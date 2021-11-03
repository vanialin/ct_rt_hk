#### responses to comments ####

require(mgcv)
##
## reviewer 2 comment 2: impact of changes in case proportion in Ct ##
## 1) symptomatic 
#ct.linelist <- read.csv("/Users/vanialam/OneDrive - connect.hku.hk/vanialam/research_vania/epi_wave_2021/program/publish/data/data_ct.csv")
#daily.linelist <- read.csv("/Users/vanialam/OneDrive - connect.hku.hk/vanialam/research_vania/epi_wave_2021/program/2021_08_R0/publish/result/daily_ct_bootstrap.csv",as.is=T)
#
data1 <- daily.linelist
# add mean age (daily) and calculate daily Ct for adult samples only
data1$mean.symp <- data1$skewness.symp <- data1$symp.record <- NA
for (i in 1:nrow(daily.linelist)){
        ct.daily <- ct.linelist[ct.linelist$date.test==data1$date[i]&
                                        ct.linelist$asymptomatic==0,]
        if (nrow(ct.daily)!=0){
                data1$mean.symp[i] <- mean(ct.daily$ct.value)
                data1$skewness.symp[i] <- 
                        e1071::skewness(ct.daily$ct.value)
        }
        data1$symp.record[i] <- nrow(ct.daily)
}
# assign data to training/testing sets
data1$period <- ifelse(as.Date(data1$date)>=as.Date("2020-07-01")&
                               as.Date(data1$date)<=as.Date("2020-08-31"),1,
                       ifelse(as.Date(data1$date)>=as.Date("2020-11-01")&
                                      as.Date(data1$date)<=as.Date("2021-03-31"),2,0))
table(data1$period) # checked
#
t.test(data1$mean[data1$period==1],data1$mean.symp[data1$period==1])
t.test(data1$skewness[data1$period==1],data1$skewness.symp[data1$period==1])
#
t.test(data1$mean[data1$period==2],data1$mean.symp[data1$period==2])
t.test(data1$skewness[data1$period==2],data1$skewness.symp[data1$period==2])
##
## get regression models (all/symptomatic) --
train.period <- seq(as.Date("2020-07-06"),as.Date("2020-08-31"),1)
train.new <- data1[data1$date%in%as.character(train.period),]
#
lm.main <- lm(log(local.rt.mean)~mean+skewness.imputed,data=train.new)
lm.symp <- lm(log(local.rt.mean)~mean.symp+skewness.symp,data=train.new)
summary(lm.symp) 
#
## get estimates for the testing period --
date.start <- as.Date("2020-11-20")
test.period <- seq(date.start,as.Date("2021-03-26"),1)
test.new <- data1[data1$date%in%as.character(test.period),]
# Ct-based Rt estimates after adjusted for daily mean age
pred.test <- exp(predict(lm.symp,test.new,interval = "prediction"))
test.new <- cbind(test.new,pred.test)
# Ct-based Rt estimates from all Ct records
pred.original <- exp(predict(lm.main,test.new,interval = "prediction"))
test.new[,c("rt.est","rt.lb","rt.ub")] <- pred.original
#  (time indicator) for plotting
test.new$test.to.start <- as.numeric(as.Date(test.new$date)-date.start)
#
#----------
## prepare plotting elements
date.seq <- seq(date.start,as.Date("2021-03-31"),1)
x.length <- length(date.seq)
month.end <- c("2020-11-30","2020-12-31","2021-01-31","2021-02-28","2021-03-31")
pos.weekend <- which(date.seq%in%c(as.Date("2020-11-20")+0:18*7))-1 # location for ticks (week)
week.end <- day(c(as.Date("2020-11-20")+0:18*7)) # text for labelling end of week
x.month <- c(3,c(which(day(date.seq)==15)-2)) # location for labelling month
x.month.pos <- c(0,(which(as.character(date.seq)%in%month.end)-1)) # for axis tick (month)
x.month.lab <- c("Nov","Dec","Jan","Feb","Mar")
#
#
#### plot out (7*12)
## panel a-b
par(mar=c(3,3,2,1)+0.1)
fig.list <- list(c(0,0.25,0.6,1),
                 c(0.5,0.75,0.6,1),
                 c(0.25,0.5,0.6,1),
                 c(0.75,1,0.6,1))
par(fig=c(0,1,0,1))
plot(NA,xlim=c(0,10),ylim=c(0,10),axes=F,type="n")
text.add <- c("a","b")
for (i in 1:2){
        df.tmp <- data1[data1$period==i,]
        par(fig=fig.list[[2*(i-1)+1]],new=T)
        boxplot(df.tmp$mean~df.tmp$rt.cat,ylim=rev(c(15,35)),axes=F,
                whisklty = 1,outpch=16,outcex=.7,staplecol="white",
                ylab=NA,xlab=NA,boxwex=.15,at=1:4-0.1,col="#ffc425")
        boxplot(df.tmp$mean.symp~df.tmp$rt.cat,ylim=rev(c(15,35)),axes=F,
                whisklty = 1,outpch=16,outcex=.7,staplecol="white",
                ylab=NA,xlab=NA,boxwex=.15,at=1:4+0.1,add=T)
        axis(1,at=1:4,labels = c(expression(""<="0.5"),"0.5-1","1-1.5",
                                 expression("">"1.5")))
        mtext("Incidence-based Rt",side=1,line=2.1)
        axis(2,at=3:7*5,las=1,line=0)
        if (i == 1){
                mtext("Daily mean Ct",side=2,line=2.1)
                mtext(text.add[1],side=3,cex=1.3,font=2,line=.5,at=0)
                polygon(c(1.1,1.1,1.3,1.3),c(15,16,16,15),col="#ffc425")
                text(1.35,15.5,"All cases",adj=0)
                polygon(c(1.1,1.1,1.3,1.3),c(17,18,18,17),col="white")
                text(1.35,17.5,"Symptomatic cases",adj=0)  
        }
        mtext(paste0("Wave ",i+2),side=3,font=2)
        # adult
        par(fig=fig.list[[2*i]],new=T)
        boxplot(df.tmp$skewness~df.tmp$rt.cat,ylim=c(-1,1.5),axes=F,
                whisklty = 1,outpch=16,outcex=.7,staplecol="white",
                ylab=NA,xlab=NA,boxwex=.15,at=1:4-0.1,
                col="#00b159")
        boxplot(df.tmp$skewness.symp~df.tmp$rt.cat,ylim=c(-1,1.5),axes=F,
                whisklty = 1,outpch=16,outcex=.7,staplecol="white",
                ylab=NA,xlab=NA,boxwex=.15,at=1:4+0.1,add=T)
        axis(1,at=1:4,labels = c(expression(""<="0.5"),"0.5-1","1-1.5",
                                 expression("">"1.5")))
        mtext("Incidence-based Rt",side=1,line=2.1)
        mtext(paste0("Wave ",i+2),side=3,font=2)
        axis(2,las=1,line=0)
        if (i == 1){
                mtext("Daily Ct skewness",side=2,line=2.5)
                mtext(text.add[2],side=3,cex=1.3,font=2,line=.5,at=0)  
                polygon(c(1.1,1.1,1.3,1.3),c(1.4,1.5,1.5,1.4),col="#00b159")
                text(1.35,1.45,"All cases",adj=0)
                polygon(c(1.1,1.1,1.3,1.3),c(1.2,1.3,1.3,1.2),col="white")
                text(1.35,1.25,"Symptomatic cases",adj=0)
        } 
}
#
# panel c
par(fig=c(0,1,0,0.6),mar=c(4,3,2,1)+0.1,new=T)
plot(NA,xlim=c(1,x.length),ylim=c(0,5),axes=F,xlab=NA)
# x-axis
day.axis <- 1:(x.length-1)
axis(1,at=day.axis,labels = rep(NA,length(day.axis)),tck=-.007)
axis(1,at=pos.weekend,labels = week.end,las=1,tck=-.012)
axis(1,at=x.month.pos,labels = rep(NA,length(x.month.pos)),tck=-.025)
for (k in 1:length(x.month)){
        mtext(x.month.lab[k],side=1,line=2.5,at=x.month[k],adj=0)
} 
# add axis legend separately
for (i in 1:2){
        axis(1,at=x.month.pos[[i]],
             labels = rep(NA,length(x.month.pos[[i]])),tck=-.06)
}
lines(test.new$test.to.start,test.new$rt.est,type="l",col="pink")
# y-axis
axis(2,at=0:5,las=1,line=-.4)
mtext("Rt",side=2,line=1.4)
lines(c(1,x.length),rep(1,2),lty=2)
polygon(c(test.new$test.to.start,
          rev(test.new$test.to.start)),
        c(test.new$rt.lb,rev(test.new$rt.ub)),
        col=alpha("pink",.3),border=F) 
# symptomatic Rt-est
for(i in 1:nrow(test.new)){
        lines(rep(test.new$test.to.start[i],2),c(test.new$upr[i],test.new$lwr[i]),
              col=alpha("#00b8ff",.45),lwd=1.7)
}
lines(test.new$test.to.start,test.new$fit,col="#00b8ff",type="p",cex=.8,pch=16)
# legends
lines(c(40,43),rep(5,2),col="pink",lwd=2)
points(41.5,5,col="pink",pch=16)
lines(c(40,43),rep(4.5,2),col="#00b8ff",lwd=2)
points(41.5,4.5,col="#00b8ff",pch=16)
text(43.5,5,"Ct predicted Rt, using all records",adj=0)
text(43.5,4.5,"Ct predicted Rt, using only symptomatic records",adj=0)
mtext("c",side=3,adj=0,font=2,cex=1.3,line=.5)
##
#####

## end of script

#####