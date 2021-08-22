#-----------
# training period
# length and timing
# Supplementary Fig. 5
#-----------
#
######################################################
## data_daily_all: daily case counts/sample counts, incidence-based Rt; 
##                 daily Ct mean, median and skewness (imputed)
##                 correspond to "Supplementary" data in source data file
######################################################
# read in "data_daily_all.csv"
daily.linelist <- read.csv("data_daily_all.csv",as.is=T)
#
# set up time period for selection training periods
startdate <- c("2020-07-01","2020-11-10")
enddate.for.plot <- as.Date(startdate)+60+20-1
#
days.select <- c(30,40,50,60) # covering XX days
r.sqr.mat <- matrix(NA,20,8)  
case.period <- list()
for (i in 1:2){ # for wave 3 and wave 4
        date.seq <- seq(as.Date(startdate[i]),as.Date(enddate.for.plot[i]),1)
        df.take <- daily.linelist[daily.linelist$date%in%as.character(date.seq),]
        case.period[[i]] <- df.take
        for (n in 1:4){ # 4 length of days covered in the training set
                for (k in 1:20){ # 20 sets each
                        df.tmp <- df.take[k:(k+days.select[n]),]
                        lm.tmp <- 
                                lm(log(local.rt.mean)~mean+skewness.imputed,
                                   data=df.tmp)
                        r.sqr.mat[k,4*(i-1)+n] <- summary(lm.tmp)$adj.r.square
                }
        }
}
#
#----------
## plot out
pdf("ED_Fig_5.pdf",height = 7,width = 11)
par(mar=c(4,4,2,3)+0.1)
fig.list <- list(c(0,0.35,0.45,1),
                 c(0,0.35,0,0.55),
                 c(0.3,1,0.45,1),
                 c(0.3,1,0,0.55))
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
dev.off()
##
#####

## end of script

#####
