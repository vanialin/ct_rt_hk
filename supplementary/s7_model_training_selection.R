#-----------
# training period
# length and timing
#-----------
#
######################################################
## data_daily_all: daily case counts/sample counts, incidence-based Rt; 
##                 daily Ct mean, median and skewness (imputed)
######################################################
require(lubridate)
require(plotrix)
# read in "data_daily_all.csv"
#daily.linelist <- read.csv("/Users/vanialam/OneDrive - connect.hku.hk/vanialam/research_vania/epi_wave_2021/program/publish/data/data_daily_all.csv",as.is=T)
#
# set up time period for selection training periods
startdate <- c("2020-07-04","2020-11-10")
enddate.for.plot <- as.Date(startdate)+60+20-1
#
days.select <- c(30,40,50,60) # covering XX days
r.sqr.mat <- case.count.mat <- matrix(NA,20,8)  
case.period <- list()
for (i in 1:2){ # for wave 3 and wave 4
        date.seq <- seq(as.Date(startdate[i]),as.Date(enddate.for.plot[i]),1)
        df.take <- daily.linelist[daily.linelist$date%in%as.character(date.seq),]
        case.period[[i]] <- df.take
        for (n in 1:4){ # 4 length of days covered in the training set
                for (k in 1:20){ # 30 sets each
                        df.tmp <- df.take[k:(k+days.select[n]),]
                        lm.tmp <- 
                                lm(log(local.rt.mean)~mean+skewness.imputed,
                                   data=df.tmp)
                        r.sqr.mat[k,4*(i-1)+n] <- summary(lm.tmp)$adj.r.square
                        case.count.mat[k,4*(i-1)+n] <- 
                                sum(df.tmp$records)/nrow(df.tmp)
                }
        }
}
#

#### added analyses
### 1. average case count
case.count.summary <- NULL
for (i in 1:nrow(r.sqr.mat)){
        for (k in 1:ncol(r.sqr.mat)){
                if (r.sqr.mat[i,k]>=0.7){
                        case.count.summary <- 
                                c(case.count.summary,case.count.mat[i,k])
                }
        }
}
summary(case.count.summary) # median (IQR) = 85 (78, 90)

### 2. using former to predict latter in wave 4
r.sqr.training <- r.sqr.mat[,5:8]
max(r.sqr.training)
which(r.sqr.training==max(r.sqr.training),arr.ind = T)

date.tmp <- seq(as.Date(startdate[2])+13,as.Date(startdate[2])+13+20-1,1)
test.tmp <- daily.linelist[daily.linelist$date%in%as.character(date.tmp),]
lm.test <- lm(log(local.rt.mean)~mean+skewness.imputed,data=test.tmp)
summary(lm.test)

## compare with main training model
train.period <- seq(as.Date("2020-07-06"),as.Date("2020-08-31"),1)
train.tmp <- daily.linelist[daily.linelist$date%in%as.character(train.period),]
lm.train <- lm(log(local.rt.mean)~mean+skewness.imputed,data=train.tmp)

####
## test with following 51 days (already cover the January peak)
pred.tmp <- daily.linelist[daily.linelist$date%in%as.character(max(date.tmp)+1:51),]
pred.tmp[,c("fit1","lwr1","upr1")] <- 
        exp(predict(lm.train,pred.tmp,interval = "prediction"))
pred.tmp[,c("fit2","lwr2","upr2")] <- 
        exp(predict(lm.test,pred.tmp,interval = "prediction"))
sum(((pred.tmp$local.rt.mean-1)*(pred.tmp$fit1-1))>=0)/nrow(pred.tmp) # 75%
sum(((pred.tmp$local.rt.mean-1)*(pred.tmp$fit2-1))>=0)/nrow(pred.tmp) # 78%

pred.tmp$date <- as.Date(pred.tmp$date)
test.tmp$date <- as.Date(test.tmp$date)
plot.tmp <- bind_rows(test.tmp,pred.tmp) %>% 
        mutate(period=1+1*(date>as.Date(max(test.tmp$date))))

#----------
## plot out
pdf("Fig_S9.pdf",height = 10,width = 10)
par(mar=c(4,4,2,3)+0.1)
fig.list <- list(c(0,0.35,0.65,1),
                 c(0,0.35,0.35,0.7),
                 c(0.3,1,0.65,1),
                 c(0.3,1,0.35,0.7))
#start.background <- c(which.max(r.sqr.mat[,1]),which.max(r.sqr.mat[,5]))
col.plot <- c("#493267","#9e379f","#e86af0","#7bb3ff")
leg.here <- c("a","b","c","d")
#
for (i in 1:2){
        if (i == 1){
                par(fig=fig.list[[i]])
        } else {
                par(fig=fig.list[[i]],new=T)
        }
        plot(NA,xlim=c(1,20),ylim=c(0,1),xlab=NA,ylab=NA,axes = F,main=NA)
        for (n in 1:4){
                lines(1:20,r.sqr.mat[,4*(i-1)+n],col=col.plot[n],lwd=2)  
                if (i == 1){
                        lines(c(5,6)+4*(n-1),rep(0.95,2),col=col.plot[n],lwd=2)
                        text(6.5+4*(n-1),0.95,days.select[n],adj=0)
                        #lines(c(16,17.5),rep(1.08-0.08*n,2),col=col.plot[n],lwd=2)
                        #text(18,1.08-0.08*n,days.select[n],adj=0)
                }
        }
        #
        date.extract <- as.Date(startdate[i])+0:19
        date.write <- date.extract[1]
        for (d in 1:20){
                if (d%%10==0){
                        date.write <- c(date.write,date.extract[d])  
                }
        }
        axis(1,at=1:20,labels = rep(NA,20),tck=-.01)
        axis(1,at=which(date.extract%in%date.write),
             labels = rep(NA,length(date.write)),tck=-.03)
        for (dd in 1:length(date.write)){
                day.tmp <- day(date.write[dd])
                month.tmp <- month(date.write[dd])
                mtext(paste0(day.tmp,"/",month.tmp),side = 1,
                      at=which(date.extract==date.write[dd])-.5,line=.5)
        }
        if (i == 2){
                mtext("Start date of the training period",side=1,line=1.7)  
        }
        #
        axis(2,las=1)
        mtext("Adjusted R square",side=2,line=2.5)
        lines(c(0.5,20.5),rep(0.7,2),lty=2,col="dark grey")
        #
        mtext(leg.here[i],side=3,adj=0,line=0,font=2,cex=1.3)
}
#
for (i in 1:2){
        df.tmp <- case.period[[i]]
        par(fig=fig.list[[i+2]],new=T)
        plot(NA,xlim=c(0,(nrow(df.tmp)-1)),ylim=c(0,150),xlab=NA,ylab=NA,
             axes = F,main=NA)
        # x-axis
        axis(1,at=0:(nrow(df.tmp)-1),tck=-0.01,
             labels = rep(NA,length(0:(nrow(df.tmp)-1))))
        date.vec <- as.character(df.tmp$date[1])
        for (d in 1:nrow(df.tmp)){
                if (d%%7==0){
                        date.vec <- c(date.vec,as.character(df.tmp$date[d]))
                }
        }
        axis(1,at=which(as.character(df.tmp$date)%in%date.vec)-1,
             labels = rep(NA,length(date.vec)),tck=-0.03)
        for (n in 1:length(date.vec)){ # x-label 
                day.tmp <- day(date.vec[n])
                month.tmp <- month(date.vec[n])
                mtext(paste0(day.tmp,"/",month.tmp),side = 1,
                      at=which(as.character(df.tmp$date)==date.vec[n])-1,line=.7)
        }
        if (i == 2){
                mtext("Date",side=1,line=1.7)
        }
        # background
        for (nn in 1:4){ # 4 alternative lengths
                start.background <- which.max(r.sqr.mat[,4*(i-1)+nn])
                polygon(c(rep(start.background-1,2),
                          rep(start.background+days.select[nn]-1,2)),
                        c(0,153-3*nn,153-3*nn,0),col=alpha(col.plot[nn],.15),border=F)
                                #
        }
        # case count and record count
        for (j in 1:nrow(df.tmp)){
                # all cases
                polygon(c(rep(j-1.5,2),rep(j-0.5,2)),c(0,rep(df.tmp$all.cases[j],2),0),
                        col="#dad5d4",border="white")
                # all records
                polygon(c(rep(j-1.5,2),rep(j-0.5,2)),c(0,rep(df.tmp$records[j],2),0),
                        col=alpha("orange",.35),border="white")
                
        }
        #
        # legend
        if (i == 1){
                polygon(c(53,53,55,55),c(150,144,144,150),
                        col="#dad5d4",border = "white")
                text(55.5,147,"Confirmed cases",adj=0)
                polygon(c(53,53,55,55),c(140,134,134,140),
                        col=alpha("orange",.3),border = "white")
                text(55.5,137,"Tested samples",adj=0)
        }
        #
        axis(2,at=0:5*30,las=1)
        mtext("Numbers",side=2,line=2.5)
        mtext(leg.here[i+2],side=3,adj=0,line=0,font=2,cex=1.3)
}
##

##

####
par(mar=c(3,4,2,4)+0.1)
par(fig=c(0,1,0,0.38),new=T)
plot(NA,xlim=c(1,nrow(plot.tmp)),ylim=c(0,150),las=1,axes=F,xlab=NA,ylab=NA)
axis(4,at=0:5*30,las=1,line=.5)
mtext("Number of cases",side=4,line=2.4)
vec <- c(0:7*10+1)
axis(1,at=vec,labels = rep(NA,5),tck=-0.02)
for (j in 1:length(vec)){
        mtext(paste0(day(as.Date(plot.tmp$date[vec[j]])),"/",
                     month(as.Date(plot.tmp$date[vec[j]]))),
              side=1,at=vec[j],line=.5) 
}
mtext("Date",side=1,line=1.2)
## plotting elements
for (i in 1:nrow(plot.tmp)){
        polygon(c(rep(i-0.5,2),rep(i+0.5,2)),
                c(0,rep(plot.tmp$records[i],2),0),col=alpha("grey",.5),border = "white")
}
par(new=T)
plot(NA,xlim=c(1,nrow(plot.tmp)),ylim=c(0,4),axes=F,xlab=NA,ylab=NA)
brackets(1,3.8,21,3.8,type = 4,xpd=T,h=.3)
mtext("Alternative training period",side=3,at=10,line=0,font=2,cex=.9)
axis(2,at=0:4,las=1,line=.5)
mtext("Rt",side=2,line=2.2)
abline(h=1,lty=2,col="grey",lwd=1.5)
polygon(c(1:nrow(plot.tmp),rev(1:nrow(plot.tmp))),
        c(plot.tmp$local.rt.lower,rev(plot.tmp$local.rt.upper)),
        col = alpha("grey",.5),border=F)
lines(1:nrow(plot.tmp),plot.tmp$local.rt.mean,lwd=2)
###
polygon(c(1:nrow(plot.tmp),rev(1:nrow(plot.tmp))),
        c(plot.tmp$lwr1,rev(plot.tmp$upr1)),
        col = alpha("pink",.3),border=F)
lines(1:nrow(plot.tmp),plot.tmp$fit1,lwd=2,col="pink")
for (i in 1:nrow(plot.tmp)){
        lines(rep(i,2),c(plot.tmp$lwr2[i],plot.tmp$upr2[i]),
              col="#7bb3ff",lwd=1.5)
        points(i,plot.tmp$fit2[i],col="#7bb3ff",pch=19,cex=.5)
}
mtext("e",side=3,font=2,adj=0,line=.5,cex=1.3)
###
## legend
col.vec <- c(1,"pink","#7bb3ff")
text.vec <- c("Incidence-based Rt",
              "Ct-based Rt, using main model",
              "Ct-based Rt, using alternative model")
for (i in 1:3){
        lines(c(35,38),rep(4.2-0.35*i,2),col=col.vec[i])
        text(38.5,4.2-0.35*i,text.vec[i],adj=0)
}

dev.off()


#####

## end of script

#####