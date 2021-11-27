options(mc.cores = 4)
require(lubridate)
require(gridExtra)
require(chron)
#source(paste0(path,"/sim_funcs_all.R"))

###
### delay distributions
set.seed(1)
incubation_period <- 
        bootstrapped_dist_fit(rlnorm(1000,log(5.2),log(3.9)),
                              dist = "lognormal",max_value = 25)
# GT for simulated data set ***
SI_mean <- pars["infectious"] + pars["incubation"]
SI_sd <- sqrt(pars["infectious"]^2 + pars["incubation"]^2)
generation_time <- list(mean=SI_mean,mean_sd=3,sd=SI_sd,sd_sd=3,max=100)

reporting_delay <- 
        bootstrapped_dist_fit(rgamma(1000,shape=1.83,rate=0.43),max_value = 25)


##
#### RUN IN WINDOWS ####
for (i in 1:4){
        ct_tmp <- read.csv(paste0(path_observe,"vl_obs1_scenario",i,".csv")) 
        get_inc_Rt(ct_tmp,n1=1,n2=i)
}

#
#

##
#### RUN IN MAC ####
### read in Rt and get results for each linelist of cases detected
load(file=paste0(path_simulate,"SEIR_dynamics",1,".Rda"))

ct_list <- rt_list <- result_list <- NULL
fit_training_list <- fit_testing_list <- box_list <- NULL
#input_range_all <- c(rep(as.Date("2020-10-19")+c(0,50),3),c(as.Date("2020-11-14")+c(0,40)))
for (i in 1:4){
        ct_tmp <- read_csv(paste0(path_observe,"vl_obs1_scenario",i,".csv"))
        rt_tmp <- read_csv(paste0(path_rt,"rt_obs1_scenario",i,".csv"))
        ct_list[[i]] <- ct_tmp
        rt_list[[i]] <- rt_tmp
        result_list[[i]] <- evaluate_daily_funcs(rt_tmp,ct_tmp,seir_dynamics)
        #
        fit_training_list[[i]] <- fit_plot(result_list[[i]],"training",add_legend = (i==1),
                                           panel = letters[i])
        fit_testing_list[[i]] <- fit_plot(result_list[[i]],"testing",
                                          input_range = c(as.Date("2020-10-19"),
                                                          as.Date("2021-01-17")))
                                                          #input_range_all[2*i]))
}

#
#

#### combine plots
##### supplementary (14*8,sup1) ----
#### plot simulation truth
## 1. case counts
case_truth <- complete_linelist%>%
        group_by(infection_time) %>%
        summarise(infect_n=n(),
                  infect_symp=sum(is_symp==1)) %>% ungroup() %>%
        ### there's no sampled time and therefore cannot be grouped by reporting dates
        ### as in this linelist the time of detection for each detected cases is not known yet
        rename(date=infection_time) %>%
        right_join(seir_dynamics$seir_outputs%>%
                           dplyr::select(step,Rt) %>% rename(date=step)) %>% 
        mutate(date=as.Date("2020-07-01")+date)  %>% arrange(date)

## 2. function for plotting dates axes
axis_plot <- function(dates_range){
        dates_all <- seq(dates_range[1],dates_range[2],1)
        dates_notation <- c(dates_all[which(days(dates_all)==1)],dates_range[2])
        axis(1,at=1:28*7,labels = rep(NA,28),tck=-.01)
        axis(1,at=which(dates_all%in%dates_notation),
             labels = rep(NA,length(dates_notation)),las=1,tck=-.025)
        for (i in 1:length(dates_notation)){
                mtext(paste0(days(dates_notation[i]),"/",
                             month(dates_notation[i]),"/",
                             substr(year(dates_notation[i]),3,4)),
                      side=1,
                      at=which(dates_all==dates_notation[i])-7,adj=0,line=.4)
        }
}

#
#

#### plot actual
## 1. count for each scenario
# a) case counts
count_list <- NULL
for (i in 1:4){
        count_list[[i]] <- 
                ct_list[[i]] %>% filter(ct_value < 40) %>% 
                group_by(sampled_time) %>% 
                summarise(count=n()) %>%
                ungroup() %>%
                right_join(ct_list[[i]]%>%group_by(infection_time) %>% 
                                   summarise(infect_count=n()) %>% ungroup() %>% 
                                   rename(sampled_time=infection_time)) %>%
                mutate(date=as.Date("2020-07-01")+sampled_time) %>%
                mutate(count=ifelse(is.na(count),0,count),
                       infect_count=ifelse(is.na(infect_count),0,infect_count))%>%
                arrange(date)
}

# b) detection proportion
detect_prob_list <- NULL
detect_prob_list[[1]] <- rep(0.25,201)
detect_prob_list[[2]] <- rep(0.1,201)
detect_prob_list[[3]] <- prob_varying$prob[1:201]
detect_prob_list[[4]] <- prob_ud$prob[1:201]

col_low <- colorRampPalette(c("#fff4e6","#be9b7b"))
col_high <- colorRampPalette(c("#be9b7b","#854442"))
col_all <- c(col_low(25),col_high(35))
#col_code <- colorRampPalette(c("#ffcc33","#ff0000"))
#col_all <- col_code(60)

#### merge all together
pdf("Fig_S9.pdf",height = 15,width = 10)
# panel a
par(fig=c(0,1,0.8,1),mar=c(3,4,1,4)+0.1)
dates_range <- range(case_truth$date)
plot(NA,xlim=c(1,nrow(case_truth)),
     ylim=c(0,round(max(case_truth$infect_n,na.rm = T),digits = -2)),
     axes=F,xlab=NA,ylab=NA)
axis_plot(dates_range = dates_range)
axis(2,las=1)
mtext("Number of cases",side=2,line=3)
mtext("Rt as simulation truth",side=4,line=2.5)
mtext("Date by infection",side=1,line=1.3)
for (i in 1:nrow(case_truth)){
        polygon(c(rep(i-0.5,2),rep(i+0.5,2)), # all infected cases
                c(0,rep(case_truth$infect_n[i],2),0),col=alpha("#eea990",.8),border="white")
        polygon(c(rep(i-0.5,2),rep(i+0.5,2)), # all symptomatic cases
                c(0,rep(case_truth$infect_symp[i],2),0),col=alpha("#aa6f73",.6),border="white")
}
### add legends
polygon(c(80,80,84,84),c(4800,5000,5000,4800),col=alpha("#eea990",.8),border="white")
polygon(c(80,80,84,84),c(4400,4600,4600,4400),col=alpha("#aa6f73",.6),border="white")
lines(c(80,84),rep(4100,2))
text(84.5,4900,"All cases",adj=0)
text(84.5,4500,"Symptomatic cases",adj=0)
text(84.5,4100,"Simulated Rt",adj=0)
##
par(new=T)
plot(NA,xlim=c(1,nrow(case_truth)),ylim=c(0,3),axes=F,xlab=NA,ylab=NA)
axis(4,las=1)
lines(c(1,nrow(case_truth)),rep(1,2),col="grey",lty=2)
lines(1:nrow(case_truth),case_truth$Rt)
mtext("a",side=3,adj=0,cex=1.2,font=2,line=-1)
#
#
## plot the color indicator
par(fig=c(0.6,1,0.74,0.8),new=T,mar=c(2,2,2,1)+0.1)
plot(NA,xlim=c(0,.6),ylim=c(0.75,1.25),axes=F,xlab=NA,ylab=NA)
axis(1,at=c(0.15,0.25,0.35,0.45,0.55),line=-1,col="white",las=1,adj=0)
#
x <- seq(0,.6,len = length(col_all))
segments(x,rep(0.7,length(x)),x,rep(1.3,length(x)), col=col_all, lwd=5.5)
mtext("Detection probability",side=3,line=.5,adj=0)
#
#
fig_list <- list(c(0,1,0.57,0.78),
                 c(0,1,0.38,0.59),
                 c(0,1,0.19,0.4),
                 c(0,1,0,0.21))
for (k in 1:4){
        #par(fig=c(0,1,0.73-0.2*(k-1),0.8-0.2*(k-1)),mar=c(2,4,2,2)+0.1,new=T)
        #plot(NA,xlim=c(1,length(dates_all)),ylim=c(0.25,0.5),axes=F,xlab=NA,ylab=NA)
        par(fig=fig_list[[k]],mar=c(3,4,2,2)+0.1,new=T)
        y_upper <- round(max(count_list[[k]]$infect_count,na.rm = T),digits = -2)
        plot(NA,xlim=c(1,nrow(case_truth)),ylim=c(0,y_upper*1.03),axes=F,xlab=NA,ylab=NA)
        axis_plot(dates_range = dates_range)
        axis(2,las=1)
        if (k == 4){
                mtext("Date",side=1,line=1.5)
        }
        if (k == 1){
                mtext("b",side=3,adj=0,cex=1.2,font=2,line=2)
                polygon(c(80,80,84,84),c(760,800,800,760),col=alpha("#aa6f73",.8),border = "white")
                text(84.5,780,"Cases by date of infection",adj=0)
                polygon(c(80,80,84,84),c(690,730,730,690),col=alpha("#bfbbbb",.85),border = "white")
                text(84.5,710,"Cases by date of reporting",adj=0)
        }
        ## plot detection
        for (i in 1:length(detect_prob_list[[k]])){
                polygon(c(rep(i-0.5,2),rep(i+0.5,2)),
                        c(y_upper,rep(y_upper*1.03,2),y_upper),
                        col=col_all[detect_prob_list[[k]][i]*100],border=NA)
        }
        
        ## plot cases
        for (i in 1:nrow(count_list[[k]])){
                polygon(c(rep(count_list[[k]]$sampled_time[i]-0.5,2),
                          rep(count_list[[k]]$sampled_time[i]+0.5,2)), 
                        # case by infection
                        c(0,rep(count_list[[k]]$infect_count[i],2),0),
                        col=alpha("#aa6f73",.8),border="white")
                polygon(c(rep(count_list[[k]]$sampled_time[i]-0.5,2),
                          rep(count_list[[k]]$sampled_time[i]+0.5,2)),
                        # case by reporting
                        c(0,rep(count_list[[k]]$count[i],2),0),
                        col=alpha("#bfbbbb",.85),border="white")
        }
        mtext("Number of cases",side=2,line=3)
        mtext(paste("scenario",k),side=3,adj = 0,line=0)
}
dev.off()


#
#

#### for output fitting results ----
mat_out <- matrix(NA,4,9)
for (i in 1:4){
        mat_out[i,1:3] <- 
                Rt_consistency(result_list[[i]]%>%
                                       filter(period==1),
                               "rt_est","fit")[2,c(1,2,4)]
        mat_out[i,4:6] <- 
                Rt_consistency(result_list[[i]]%>%
                                       filter(period==2),
                               "rt_est","fit")[2,c(1,2,4)]
        mat_out[i,7:9] <- 
                Rt_consistency(result_list[[i]]%>%
                                       filter(period==2&
                                                      count>30),
                               "rt_est","fit")[2,c(1,2,4)]
}

#write.csv(mat_out,"table_s5.csv",row.names = F)
rho_vec <- mat_out[1:4,4]
for (i in 1:4){
        box_list[[i]] <- boxplot_func(result_list[[i]]%>%filter(period==2),
                                      add_legend = i==1,
                                      rho_ct = rho_vec[i])
}


# results (14*20,result1)
grid.arrange(grobs=list(
        fit_training_list[[1]],fit_testing_list[[1]],box_list[[1]],
        fit_training_list[[2]],fit_testing_list[[2]],box_list[[2]],
        fit_training_list[[3]],fit_testing_list[[3]],box_list[[3]],
        fit_training_list[[4]],fit_testing_list[[4]],box_list[[4]]), ## example of how to arrange it
        heights = rep(1,4),
        widths = c(2,4,2),
        layout_matrix = rbind(c(1,2,3),
                              c(4,5,6),
                              c(7,8,9),
                              c(10,11,12)))
        
Rt_consistency(result_list[[1]]%>%filter(period==1),"rt_est","fit")[2,]
Rt_consistency(result_list[[1]]%>%filter(period==2),"rt_est","fit")[2,]




mat_tmp <- matrix(NA,8,8)
for (i in 1:4){
        mat_tmp[2*(i-1)+1:2,1:4] <- 
                Rt_consistency(result_list[[i]]%>%
                                       filter(period==2),
                                       "rt_est","fit")
        mat_tmp[2*(i-1)+1:2,5:8] <- 
                Rt_consistency(result_list[[i]]%>%
                                      filter(period==2&count>30),
                               "rt_est","fit")
}





####
dates_check <- seq(as.Date("2020-10-19"),as.Date("2020-11-14"),1)
dates_check <- as.Date("2020-07-01")+60:110
boxplot_func(result_list[[4]]%>%filter(date>as.Date("2020-10-19")),T)

Rt_consistency(result_list[[4]]%>%filter(date%in%dates_check),"rt_est","fit")



