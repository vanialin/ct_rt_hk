#------------
# back-projecting Ct
# on illness onset
#------------
#
# load packages
require(e1071) # to calculate skewness
require(lubridate)
#
######################################################
## data_ct: all individual Ct values (with test dates)
## data_cases: daily case counts/sample counts and incidence-based Rt
######################################################
# read in "data_ct.csv" and "data_cases.csv"
#ct.linelist <- read.csv("/Users/vanialam/OneDrive - connect.hku.hk/vanialam/research_vania/epi_wave_2021/program/2021_08_R0/publish/data/data_ct.csv")
#daily.linelist <- read.csv("/Users/vanialam/OneDrive - connect.hku.hk/vanialam/research_vania/epi_wave_2021/program/2021_08_R0/publish/data/data_cases.csv")
#
ct.symp <- ct.linelist[!is.na(ct.linelist$date.onset),]
# log-linear trend of Ct over time
# key - individual indicator (i.e., patient ID)
lm1 <- lm(log(ct.value)~test.to.onset*age.gp+factor(key),data=ct.symp)
coef.grab <- coef(lm1)
coef.df <- data.frame(key.factor=names(coef.grab),value=as.numeric(coef.grab))
coef.df$ct.onset <- coef.df$age.gp <- coef.df$key <- NA
for (i in 4:nrow(coef.df)){
        coef.df$key[i] <- substr(as.character(coef.df$key.factor[i]),12,16)
        coef.df$age.gp[i] <- ct.symp$age.gp[ct.symp$key==coef.df$key[i]][1]
        coef.df$ct.onset[i] <- 
                exp(coef.df$value[1]+
                            coef.df$value[3]*coef.df$age.gp[i]+
                            coef.df$value[i])
}
tail(coef.df) # check
coef.out <- coef.df[4:(nrow(coef.df)-1),]
# alternative method to calculate the Ct at onset
#coef.out$test.to.onset <- 0
#coef.out$ct.onset <- exp(predict(lm1,coef.out)) 
# 
coef.out <- coef.out[,c("key","age.gp","ct.onset")]
#
ct.linelist2 <- merge(ct.linelist,coef.out[,c("key","ct.onset")],
                      by.x = "key",all.x = T,by.y = "key")
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
#
##
#########
## show skewness and value
## between onset and actual
#########
## prepare plotting elements
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
                onset.used <-  #categorized on DATE.ONSET 
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
                    min(data.tmp$local.rt.mean)<=1){ # half days' declining
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
#---------
## start plotting
#pdf("ED_Fig_2.pdf",height = 10,width = 9)
par(mar=c(4,3,2,2)+0.1)
#x.vec <- c(c(0.25,0.5,0.75,1),rep(c(0.25,0.5,0.75),2))
x.vec <- c(rep(0.4,4),rep(0.7,3),rep(1,3))
#y.vec <- c(rep(1,4),rep(c(0.67,0.4),each=3))
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
        #if (i %in% 8:10){
         #       par(mar=c(3,3.5,4,0.5)+0.1) 
        #} else {
        #        par(mar=c(4,3.5,3,0.5)+0.1)
        #}
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
# increasing
text(0,6.5,"Decreasing epidemic",adj=0)
polygon(c(rep(0,2),rep(1,2)),c(5.5,5.7,5.7,5.5),col=col.option[1], # all record by sampling
        border="white")
text(1.3,5.6,"Ct at sampling",adj=0)
polygon(c(rep(0,2),rep(1,2)),c(4.5,4.7,4.7,4.5),col=alpha(col.option2[1],.7), # all record by onset
        border="white")
text(1.3,4.6,"Ct at onset",adj=0)
dev.off()
##
##-----
# calculate coefficients of variation
sd(skew.actual)/mean(skew.actual)
sd(skew.onset)/mean(skew.onset)
##
#####

## end of script

#####