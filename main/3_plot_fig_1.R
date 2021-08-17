#------------
# plot Ct distributions
# Fig. 1 in paper
#-------------
######################################################
## daily_ct_bootstrap: daily case counts/sample counts, incidence-based Rt; 
##                     daily Ct mean, median and skewness (imputed)
##                     CIs calculated for GAM Ct and skewness (from "2_ct_for_bootstrap")
######################################################
# load packages
require(plotrix)
require(lubridate)
#
# read in "daily_ct_bootstrap.csv"
daily.ct <- read.csv("daily_ct_bootstrap.csv")
#
#
start.date <- c("2020-07-01","2020-11-01")
end.date <- c("2020-08-31","2021-03-31")
# set up plotting elements
# x-axis
x.day.list <- x.month.pos <- x.month.list <- x.pos.list <- list()
x.length <- rep(NA,2)
month.end <- list(c("2020-07-31","2020-08-31"),
                  c("2020-11-30","2020-12-31","2021-01-31","2021-02-28","2021-03-31"))
for (i in 1:2){
        # start of each 7-day interval
        if (i == 1){
                x.tmp <- c(as.Date(start.date[i])+1:8*7)
        } else {
                x.tmp <- c(as.Date(start.date[i])+1:21*7)
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
#--------
## start plotting
pdf("Fig_1.pdf",width = 15,height = 8)
## panel A: epi curve
par(mar=c(4.5,3,3,3)+0.1)
par(fig=c(0,0.7,0.6,1))
plot(NA,xlim=c(1,x.total),ylim=c(0,150),xlab=NA,ylab=NA,axes = F,main=NA)
# y-axis for number of cases
axis(4,0:5*30,labels=rep(NA,6),line=-1.5)
mtext(0:5*30,side=4,at=0:5*30,las=1,line=-.5)
mtext("Number of cases",side=4,line=1.4)
# x-axis
day.axis <- c(0:(x.length[1]-1),gap.wave4:(x.length[2]+gap.wave4-1))
axis(1,at=day.axis,
     labels = rep(NA,length(day.axis)),tck=-.02)
axis(1,at=x.pos.combine,labels = rep(NA,length(x.pos.combine)),las=1,tck=-.04)
axis(1,at=c(0,131,221),labels=rep(NA,3),tck=-.13)
axis.break(1,x.length[1]+4,style="slash",brw=.015)
# add axis legend separately
for (i in 1:2){
        axis(1,at=x.month.pos[[i]],
             labels = rep(NA,length(x.month.pos[[i]])),tck=-.06)
}
# case count by reporting date
for (i in 1:2){
        df.tmp <- daily.ct[daily.ct$period==i,]
        for (j in 1:nrow(df.tmp)){
                polygon(c(rep(j-1.5+start.vec[i],2),rep(j-0.5+start.vec[i],2)),
                        c(0,rep(df.tmp$all.cases[j],2),0),
                        col="#dad5d4", # all cases
                        border="white")
        }
}
## Rt part
par(new=T)
plot(NA, ylim=c(0,5),xlim=c(0,x.total),xlab=NA,ylab=NA,axes=F)
axis(2,at=0:5,las=1,line=-.8)
mtext("Incidence-based Rt",side=2,line=1)
lines(c(0,x.total),rep(1,2),lty=2)
for (i in 1:2){
        rt.tmp <- daily.ct[daily.ct$period==i,]
        polygon(c(rt.tmp$test.to.start+start.vec[i],
                  rev(rt.tmp$test.to.start)+start.vec[i]),
                c(rt.tmp$local.rt.lower,rev(rt.tmp$local.rt.upper)),
                col=alpha("black",.2),border=F) 
        lines(rt.tmp$test.to.start+start.vec[i],rt.tmp$local.rt.mean,
              col=alpha("black",.9),lwd=2) # Rt (actual)
}
# add bracket
for (i in 1:2){
        brackets(start.vec[i],5,
                 start.vec[i]+x.length[i],5,type = 4,xpd=T,h=.3)
        mtext(month.title[i],side=3,at=x.length[i]/2+start.vec[i],
              line=.8,font=2,cex=1)
}
# legend
polygon(c(160,160,163,163),c(4.4,4.6,4.6,4.4),col="#dad5d4",border="white")
text(163.5,4.5,"Cases by reporting date",adj=0)
lines(c(160,163),rep(3.8,2),lwd=2)
text(163.5,3.8,"Incidence-based Rt",adj=0)
mtext("a",side=3,line=.5,adj=0,font=2,cex=1.4)
##
## panel B: Ct distribution
par(fig=c(0,0.7,0.3,0.7),new=T)
plot(NA,xlim=c(1,x.total),ylim=c(0,150),xlab=NA,ylab=NA,axes = F,main=NA)
# y-axis for number of cases
axis(4,0:5*30,labels=rep(NA,6),line=-1.5)
mtext(0:5*30,side=4,at=0:5*30,las=1,line=-.5)
mtext("Number of records",side=4,line=1)
# x-axis
day.axis <- c(0:(x.length[1]-1),gap.wave4:(x.length[2]+gap.wave4-1))
axis(1,at=day.axis,
     labels = rep(NA,length(day.axis)),tck=-.02)
axis(1,at=x.pos.combine,labels = rep(NA,length(x.pos.combine)),
     las=1,tck=-.04,cex.axis=1)
axis(1,at=c(0,131,221),labels=rep(NA,3),tck=-.13)
axis.break(1,x.length[1]+4,style="slash",brw=.015)
# add axis legend separately
for (i in 1:2){
        axis(1,at=x.month.pos[[i]],
             labels = rep(NA,length(x.month.pos[[i]])),tck=-.06)
}
# record count by sampling date
for (i in 1:2){
        df.tmp <- daily.ct[daily.ct$period==i,]
        for (j in 1:nrow(df.tmp)){
                polygon(c(rep(j-1.5+start.vec[i],2),rep(j-0.5+start.vec[i],2)),
                        c(0,rep(df.tmp$records[j],2),0),
                        col="#838584", # all record
                        border="white")
        }
}
par(new=T)
plot(NA,xlim=c(1,x.total),ylim=rev(c(18,28)),xlab=NA,ylab=NA,axes = F,main=NA)
axis(2,at=9:14*2,las=1,cex.axis=1,line=-.8)
mtext("Ct value",side=2,line=1.5)
# plot GAM
for (i in 1:2){
        df.tmp <- daily.ct[daily.ct$period==i,]
        polygon(c(df.tmp$test.to.start+start.vec[i],
                  rev(df.tmp$test.to.start)+start.vec[i]),
                c(df.tmp$gam.lb,rev(df.tmp$gam.ub)),
                col=alpha("orange",.3),border=F) 
        lines(df.tmp$test.to.start+start.vec[i],df.tmp$ct.gam,col="orange",lwd=1.5)
}
# legend
polygon(c(160,160,163,163),c(18.5,19,19,18.5),col="#838584",border="white")
text(163.5,18.75,"Records by sampling date",adj=0)
lines(c(160,163),rep(20,2),lwd=2,col="orange")
text(163.5,20,"Ct smoothed by GAM",adj=0)
mtext("b",side=3,line=.5,adj=0,font=2,cex=1.4)
##
## panel C: Ct skewness
### skewness
par(fig=c(0,0.7,0,0.4),new=T)
plot(NA,xlim=c(1,x.total),ylim=c(-2,2),xlab=NA,ylab=NA,axes = F,main=NA)
axis(2,at=-2:2,las=1,line=-.8)
mtext("Skewness",side=2,line=1.5)
lines(c(1,x.total),rep(0,2),col="grey",lty=2)
# x-axis
day.axis <- c(0:(x.length[1]-1),gap.wave4:(x.length[2]+gap.wave4-1))
axis(1,at=day.axis,
     labels = rep(NA,length(day.axis)),tck=-.02)
axis(1,at=x.pos.combine,labels = x.day.combine,las=1,tck=-.04)
axis(1,at=c(0,131,221),labels=rep(NA,3),tck=-.13)
axis.break(1,x.length[1]+4,style="slash",brw=.015)
# add axis legend separately
for (i in 1:2){
        axis(1,at=x.month.pos[[i]],
             labels = rep(NA,length(x.month.pos[[i]])),tck=-.06)
        for (k in 1:length(x.month.list[[i]])){
                mtext(x.month.lab[[i]][k],side=1,line=2.3,
                      at=x.month.list[[i]][k],adj=0)
        } 
}
# add longer axis
mtext("2020",at = 66, side=1,line=3.5,font=2)
mtext("2021",at = 175, side=1,line=3.5,font=2)
# skewness
for (i in 1:2){
        df.tmp <- daily.ct[daily.ct$period==i,]
        for (n in 1:nrow(df.tmp)){
                if (!(is.na(df.tmp$skewness[n])|is.nan(df.tmp$skewness[n]))){
                        lines(rep(df.tmp$test.to.start[n]+start.vec[i],2),
                              c(df.tmp$skewness.ub[n],df.tmp$skewness.lb[n]),
                              col=alpha("#005b96",.6))
                        points(df.tmp$test.to.start[n]+start.vec[i],
                               df.tmp$skewness[n],col="#005b96",pch=18,cex=.7)
                }
        }
}
mtext("c",side=3,line=.5,adj=0,font=2,cex=1.4)
####
####
## correlation panels
## panel D: mean Ct and Rt
df.list <- list()
for (i in 1:2){
        df.list[[i]] <- daily.ct[daily.ct$period==i,]
}
par(fig=c(0.7,1,0.5,1),mar=c(5,4,4,1)+0.1,new=T)
plot(NA,xlim=c(0.5,5.5),ylim=c(0,4),xlab=NA,ylab=NA,axes=F)
lines(c(0.5,5.5),rep(1,2),lty=2,col="grey")
boxplot(df.list[[1]]$local.rt.mean~df.list[[1]]$mean.cat,ylim=c(0,4),axes=F,
        ylab=NA,xlab=NA,boxwex=.15,at=1:5-0.1,whisklty = 1,outpch=16,outcex=.7,staplecol="white",
        col="orange",add=T)
boxplot(df.list[[2]]$local.rt.mean~df.list[[2]]$mean.cat,ylim=c(0,4),axes=F,
        ylab=NA,xlab=NA,boxwex=.15,at=1:5+0.1,whisklty = 1,outpch=16,outcex=.7,staplecol="white",
        col="#ffdbac",add=T)
axis(1,at=1:5,labels = c(expression(""<="20"),"20-22","22-24","24-26",
                         expression("">"26")))
mtext("Incidence-based Rt",side=2,line=1)
mtext("Daily mean Ct",side=1,line=2.5)
axis(2,at=0:4,las=1,line=-1)
## add legend
polygon(c(3.2,3.2,3.4,3.4),c(3.7,3.8,3.8,3.7),col="orange")
text(3.45,3.75,"Training",adj=0)
polygon(c(3.2,3.2,3.4,3.4),c(3.3,3.4,3.4,3.3),col="#ffdbac")
text(3.45,3.35,"Testing",adj=0)
mtext("d",side=3,line=1,adj=0,font=2,cex=1.4)
##
## panel E: skewness and Rt
par(fig=c(0.7,1,0,.5),new=T)
plot(NA,xlim=c(0.5,4.5),ylim=c(0,4),xlab=NA,ylab=NA,axes=F)
lines(c(0.5,4.5),rep(1,2),lty=2,col="grey")
boxplot(df.list[[1]]$local.rt.mean~df.list[[1]]$skewness.cat,ylim=c(0,4),axes=F,
        ylab=NA,xlab=NA,boxwex=.15,at=1:4-0.1,whisklty = 1,outpch=16,outcex=.7,staplecol="white",
        col="#005b96",add=T)
boxplot(df.list[[2]]$local.rt.mean~df.list[[2]]$skewness.cat,ylim=c(0,4),axes=F,
        ylab=NA,xlab=NA,boxwex=.15,at=1:4+0.1,whisklty = 1,outpch=16,outcex=.7,staplecol="white",
        col="white",add=T)
axis(1,at=1:4,labels = c(expression(""<="-0.3"),"(-0.3,0]","(0,0.3]",
                         expression("">"0.3")))
mtext("Daily Ct skewness",side=1,line=2.5)
mtext("Incidence-based Rt",side=2,line=1)
axis(2,at=0:4,las=1,line=-1)
## add legend
polygon(c(2.6,2.6,2.8,2.8),c(3.7,3.8,3.8,3.7),col="#005b96")
text(2.85,3.75,"Training",adj=0)
polygon(c(2.6,2.6,2.8,2.8),c(3.3,3.4,3.4,3.3),col="white")
text(2.85,3.35,"Testing",adj=0)
mtext("e",side=3,line=1,adj=0,font=2,cex=1.4)
#
dev.off()
##
#####

## end of script

#####