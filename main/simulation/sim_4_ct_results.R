#------------
# calculate daily Ct and merge with inc-based Rt
# build regression and select training period
# output results
#------------

## load packages
require(mgcv)
require(e1071)
require(gridExtra)
require(grid)
require(DescTools) # for Spearman CI
require(ggpubr)

source(paste0(path,"/sim_funcs_others.R"))
source(paste0(path,"/sim_funcs_evaluate_simple.R"))

load(file=paste0(path,"/linelists/SEIR_dynamics.Rda"))
#simulate_truth <- seir_dynamics$seir_outputs # simulated "true" Rt

### read in data and store result for each linelist directly ----
result_list <- NULL
for (i in 1:9){
        rt_tmp <- read.csv(paste0(path_output,"rt_obs1_scenario",i,".csv"), as.is = T)
        ct_tmp <- read.csv(paste0(path_observe,"vl_obs1_scenario",i,".csv"))
        
        result_list[[i]] <- 
                evaluate_daily_funcs(rt_tmp,ct_tmp,sim_dyn=seir_dynamics, # input data
                                     ## demanded by nested function
                                     add_legend = 1*(i==1))
        
        ### OUTPUTS HERE:
        ## 1) Epi-curve and smoothing Ct (for preliminary check)
        ## 2) correlation between Ct and Rt (training vs. testing)
        ## 3) plot for selected training period (based on adjusted R square)
        ## 4) model outputs from model built upon the best training period
        ## 5) boxplot comparing consistency between estimates and simulation truth (similar to Fig. 3b)
        ## 6) Epi-curve and all Rt (inc-based, Ct-based, truth) over training and testing (similar to Fig. 3a)
        ## 7) the whole dataframe with all information
}

#
#

#### OUTPUTS ####
### output 1: accuracy of estimates (current Fig.4) (8*10) ----
par(mar=c(3,3,2,2)+0.1)
par(fig=c(0,0.5,0.5,1))
plot(NA,xlim=c(0.5,9.5),ylim=c(-1,1.2),las=1,axes=F)
axis(1,at=1:9)
axis(2,at=(-2:2)*0.5,las=1)
mtext("Coefficients",side=2,line=2.2)
mtext("Scenario",side=1,line=2)
lines(c(0.5,9.5),rep(0,2),lty=3)
col.mean <- c(1,"dark grey")
col.skew <- c("dark blue","#b3cde0")
for (i in 1:9){
        if (i %% 2 !=0){
                polygon(c(rep(i-0.5,2),rep(i+0.5,2)),
                        c(-1,1,1,-1),col=alpha("grey",.3),border = F)
        }
        df_tmp <- result_list[[i]]$df_out
        train <- df_tmp %>% filter(period==1)
        rho.mean <- SpearmanRho(log(train$rt_est),train$mean,
                                use="na.or.complete",conf.level = 0.95)
        rho.skewness <- SpearmanRho(log(train$rt_est),train$skewness,
                                    use="na.or.complete",conf.level = 0.95)
        for_ci_plot(i-0.15,rho.mean,col=1)
        for_ci_plot(i+0.15,rho.skewness,col="dark blue")
}
lines(c(6.2,6.5),rep(1.2,2),col=1,lwd=1.5);lines(rep(6.35,2),c(1.18,1.22),col=1,lwd=2)
text(6.7,1.2,"Mean",adj=0)
lines(c(6.2,6.5),rep(1.1,2),col="dark blue",lwd=1.5)
lines(rep(6.35,2),c(1.08,1.12),col="dark blue",lwd=2)
text(6.7,1.1,"Skewness",adj=0,lwd=2)
mtext("a",side=3,font=2,adj=0,cex=1.3)

##
par(fig=c(0.5,1,0.5,1),new=T)
plot(NA,xlim=c(0.5,9.5),ylim=c(0,1.1),las=1,axes=F)
axis(1,at=1:9)
axis(2,at=(0:4)*0.25,las=1)
mtext("Spearman Rho",side=2,line=2.7)
mtext("Scenario",side=1,line=2)
lines(c(0.5,9.5),rep(0.5,2),lty=3)
col.list <- list(c("#234d20","#c9df8a"),
                 c("dark red","#f7d0cb"))
for (i in 1:9){
        if (i %% 2 !=0){
                polygon(c(rep(i-0.5,2),rep(i+0.5,2)),
                        c(0,1,1,0),col=alpha("grey",.3),border = F)
        }
        df_tmp <- result_list[[i]]$df_out
        for (j in 1:2){
                rho <- round(Rt_consistency(df_tmp %>% filter(period==j),
                                             "rt_est","fit")[2,1:3],2)
                for_ci_plot(i-0.15+0.3*(j-1),rho,col=col.list[[j]][1])
        }
}
lines(c(6.2,6.5),rep(1.1,2),col="#234d20",lwd=1.5)
lines(rep(6.35,2),c(1.09,1.11),col="#234d20",lwd=2)
text(6.7,1.1,"Training",adj=0)
lines(c(6.2,6.5),rep(1.05,2),col="dark red",lwd=1.5)
lines(rep(6.35,2),c(1.04,1.06),col="dark red",lwd=2)
text(6.7,1.05,"Testing",adj=0,lwd=2)
mtext("b",side=3,adj=0,cex=1.3,font=2)

##
## based on case count 
par(fig=c(0,1,0,0.55),new=T)
col.used <- c("#4ebcff","#2972b6","#002790","#945cb4","#001d4f")
text.used <- c("<100","100-500","500-1000","1000-5000",">5000")
plot(NA,xlim=c(0,6*9),ylim=c(0,160),axes=F,xlab="Scenario",ylab=NA)
axis(1,at=0:8*6+3,labels = 1:9)
axis(4,las=1,at=0:7*20,line=-1)
mtext("Days",side=4,line=1)
mtext("Scenario",side=1,line=2)
for (i in 1:9){
        df <- result_list[[i]]$df_out %>%
                mutate(case_cut = cut(count,
                                      breaks = c(0,100,500,1000,5000,10070),
                                      labels = 1:5))
        for (k in 1:5){
                polygon(c(rep(6*(i-1)+k-0.5,2),rep(6*(i-1)+k+0.5,2)),
                        c(0,rep(nrow(df[df$case_cut==k,]),2),0),
                        col=alpha(col.used[k],.5),border="white")
                polygon(c(rep(-5+10*k+2*(k==1),2),rep(-3+10*k+2*(k==1),2)),
                        c(143,147,147,143),col=alpha(col.used[k],.3),border="white")
                text(-2.7+10*k+2*(k==1),145,text.used[k],adj=0)
        }
}
par(new=T)
plot(NA,xlim=c(0,6*9),ylim=c(-0.2,1.2),axes=F,xlab=NA,ylab=NA)
axis(2,las=1,at=-1:5*0.2,line = -1)
mtext("Spearman Rho",side=2,line=1)
lines(c(0,54),rep(0,2),col="grey",lty=2)
lines(c(0,54),rep(0.5,2),col="grey",lty=2)
for (i in 1:9){
        df <- result_list[[i]]$df_out %>%
                filter(period == 2) %>%
                mutate(case_cut = cut(count,
                                      breaks = c(0,100,500,1000,5000,10070),
                                      labels = 1:5))
        for (k in 1:5){
                if (nrow(df[df$case_cut==k,])!=0){
                        cor <- Rt_consistency(df[df$case_cut==k,],"rt_est","fit")[2,1:3]
                        for_ci_plot(6*(i-1)+k,cor,col="dark red")
                }
                
        }
}
mtext("c",side=3,font=2,adj=0,cex=1.3,line=-1.5)

#
#

### output 2: comparison of inc-Rt under all scenarios (Fig S11)----
p_list <- NULL
for (i in 1:9){
     p_list[[i]] <- 
             epi_plot(result_list[[i]]$df_out,truth = T,training = T)$epi + 
             theme(legend.position = c(0.8, 0.85),
                   axis.title = element_text(size = 10),
                   plot.title = element_text(size = 10, face = 'bold')) +
             labs(title = paste0(letters[i]," - scenario ",i)) 
}
#
(p <- grid.arrange(grobs=list(p_list[[1]],p_list[[2]],p_list[[3]],
                              p_list[[4]],p_list[[5]],p_list[[6]],
                              p_list[[7]],p_list[[8]],p_list[[9]]),
                   heights = c(1, 1, 1),
                   widths = c(1,1,1),
                   layout_matrix = rbind(c(1,2,3),
                                         c(4,5,6),
                                         c(7,8,9)))) # 9*12

#
#

### output 3: performance under symptom-based setting (Fig 3)(5*10) ----
plot_left <- rt_plot(result_list[[1]]$df_out,time="all",add_epi = T) +
        theme(plot.title = element_text(face = 'bold')) +
        labs(title = "a") 
plot_right <- result_list[[1]]$plot2 + 
        theme(plot.title = element_text(face = 'bold')) +
        labs(title = "b") 
(q <- grid.arrange(grobs=list(plot_left,plot_right),
                  heights = 1,
                  widths = c(5,2),
                  layout_matrix = rbind(c(1,2))))

#
#

### output 4: directional consistency (Table S5)
dir_mat <- matrix(NA,3,9)
tmp_list <- list(c(1,2),1,2)
for (k in 1:3){ # all/training/testing
        for (i in 1:9){
                dir_mat[k,i] <- 
                        Rt_consistency(result_list[[i]]$df_out %>% 
                                                filter(period %in% tmp_list[[k]]),
                                        "rt_est","fit")[2,5]
        }  
}

#write.csv(dir_mat,paste0(path_plot,"directional_31oct.csv"),row.names = F)

### output 5: compare performance under varying detection (table s6) ----
# reference is other scenario with flat and similar detection proportions
compare_mat <- matrix(NA,10,6)
compare_seq <- list(c(2,7,8),c(3,4,5),c(3,9))
period_list <- list(c(1,2),1,2)
for (k in 1:3){ # 3 comparison sets 
        compare_tmp <- compare_seq[[k]]
        for (i in 1:length(compare_tmp)){ # each compare group
                for (j in 1:3){ # all,training, testing
                        df_tmp <- 
                                result_list[[compare_tmp[i]]]$df_out%>%
                                filter(period%in%period_list[[j]])
                        a <- round(Rt_consistency(df_tmp,"rt_est","fit"),2)
                        compare_mat[4*(k-1)+i,2*(j-1)+1] <- paste0(a[2,1],"(",a[2,2],",",
                                                                   a[2,3],")")
                        compare_mat[4*(k-1)+i,2*j] <- a[2,5]
                }
                
        }
        
}

compare_mat[is.na(compare_mat)] <- ""
#write.csv(compare_mat,paste0(path_plot,"table_s6.csv"),row.names = F)

##
#####

## end of script

#####