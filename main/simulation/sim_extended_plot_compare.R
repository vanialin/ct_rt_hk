#-----------
# plot to check
# re-fitted viral kinetics parameters
#-----------
#
## Standard deviations of viral kinetics priors
sds_exp  <- c("beta"=0.25,"R0"=0.6,
              "obs_sd"=0.5,"viral_peak"=2,
              "wane_rate2"=1,"t_switch"=3,"level_switch"=1,
              "prob_detect"=0.03,
              "incubation"=0.25, "infectious"=0.5,
              "rho"=2,"nu"=0.5)

## ggplot sources
export_theme <- theme_tufte() +
        theme(
                axis.text.x = element_text(size=7,family="sans"),
                axis.text.y=element_text(size=7,family="sans"),
                axis.title.x=element_text(size=8,family="sans",vjust=-1),
                axis.title.y=element_text(size=8,family="sans"),
                
                ## Axis lines
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                
                ## Title
                plot.title = element_text(family="sans",size=8,face="bold",hjust=0.5),
                plot.tag = element_text(family="sans",size=10,face="bold"),
                
                ## Legends
                legend.title=element_text(size=8,family="sans",face="italic"),
                legend.text=element_text(size=8,family="sans"),
                legend.key.size= unit(0.5, "cm"),
                legend.margin = margin(0,0,0,0, "cm"),
                ## Strips for facet_wrap
                strip.text=element_text(size=8,family="sans",face="bold"),
                strip.background=element_blank()
        )


###
parTab1 <- read.csv(paste0(path_source,"partab_fitted_hk.csv"))
parTab2 <- read.csv(paste0(path_source,"partab_fitted_nh.csv"))

pars_list <- list(parTab1,parTab2)

p_main <- NULL
for (compare in 1:2){
        parTab_used <- pars_list[[compare]]
        pars_used <- parTab_used$values
        names(pars_used) <- parTab_used$names
        #
        #
        test_ages <- seq(0,35,by=1)
        vls <- viral_load_func(pars_used, test_ages, FALSE)
        
        ## Generate draws from our priors and get quantiles
        sample_pars_from_prior <- function(pars){
                sds <- sds_exp
                ## Draw parameters from normal priors
                tmp_pars <- pars
                tmp_pars["obs_sd"] <- max(rnorm(1, pars["obs_sd"], sds["obs_sd"]), 1)
                tmp_pars["viral_peak"] <- rnorm(1, pars["viral_peak"], sds["viral_peak"])
                tmp_pars["t_switch"] <- rnorm(1,pars["t_switch"],sds["t_switch"])
                
                beta1_mean <- pars["prob_detect"]
                beta1_sd <- sds["prob_detect"]
                beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
                beta_beta <- beta_alpha*(1/beta1_mean - 1)
                tmp_pars["prob_detect"] <- rbeta(1,beta_alpha,beta_beta)
                tmp_pars
        }
        
        ## Draw 10000 samples from prior
        n_samp <- 10000
        tmp_trajs <- matrix(0, nrow=n_samp, ncol=length(test_ages))
        all_pars <- matrix(nrow=n_samp,ncol=length(pars_used))
        
        for(i in 1:n_samp){
                ## Get parameters
                tmp_pars <- sample_pars_from_prior(pars_used)
                ## Get trajectory Ct scale
                tmp_vls <- viral_load_func(tmp_pars, test_ages, FALSE)
                ## Save
                tmp_trajs[i,] <- tmp_vls
                all_pars[i,] <- tmp_pars
        }
        
        ## Get median, 95 and 50% quantiles
        vl_quants <- t(apply(tmp_trajs, 2, function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975))))
        vl_quants <- as.data.frame(vl_quants)
        colnames(vl_quants) <- c("lower","mid_low","median","mid_high","upper")
        vl_quants$t <- test_ages
        
        tmp_vl_melted <- reshape2::melt(tmp_trajs)
        colnames(tmp_vl_melted) <- c("samp","t","vl")
        tmp_vl_melted$t <- test_ages[tmp_vl_melted$t]
        
        ## trajectory of shedding ##
        p_vl <- ggplot(vl_quants) + 
                geom_ribbon(aes(x=t,ymin=lower,ymax=upper), alpha=0.25,fill="#D55E00") +
                geom_ribbon(aes(x=t,ymin=mid_low,ymax=mid_high),alpha=0.5,fill="#D55E00") +
                scale_y_continuous(trans="reverse") +
                coord_cartesian(ylim=c(40,10)) +
                geom_line(aes(x=t,y=median),col="#D55E00",lwd=1.2) + 
                theme(plot.tag = element_text(size=13,face="bold")) +
                labs(tag=letters[2*(compare-1)+1])
        
        #
        ## distributions of core parameters ## 
        ## MORE samples from prior
        n_samps <- 100000
        ## Get point estimate trajectory
        vls_est <- viral_load_func(pars_used, test_ages, FALSE)
        line_dat <- tibble(mean_load=vls_est,t=test_ages)
        
        colnames(all_pars) <- parTab_used$names
        all_pars_melted <- reshape2::melt(all_pars)
        colnames(all_pars_melted) <- c("samp","par","value")
        parTab_used[parTab_used$names == "beta","fixed"] <- 1
        all_pars_melted <- all_pars_melted %>% 
                filter(par %in% parTab_used[parTab_used$fixed == 0, "names"])
        
        par_key <- c("obs_sd"="sigma[obs]",
                     "viral_peak"="VL[peak]",
                     "t_switch"="t[switch]",
                     "level_switch"="VL[switch]",
                     "wane_rate2"="t[LOD]",
                     "prob_detect"="p[addl]"
        )
        all_pars_melted$par <- as.character(all_pars_melted$par)
        all_pars_melted$label <- par_key[all_pars_melted$par]
        all_pars_melted$label <- factor(all_pars_melted$label,
                                        levels=c("sigma[obs]","VL[peak]","t[switch]",
                                                 "VL[switch]","t[LOD]","p[addl]"))
        
        p_dist <- 
                ggplot(all_pars_melted %>% filter(par %in% c("viral_peak","obs_sd","t_switch"))) +
                geom_density(aes(x=value), fill="#0072B2",alpha=0.5) +
                ylab("Prior density") +
                xlab("Value") +
                facet_wrap(~label, scales="free", ncol=3,labeller=label_parsed) +
                export_theme+
                theme(axis.text.x=element_text(size=6),
                      axis.text.y=element_text(size=6),
                      plot.tag = element_text(size=13)) +
                labs(tag=letters[2*compare])
        p_main[[compare]] <- (p_vl | p_dist) + plot_layout(widths = c(1,1.5)) 
}


p_out <- (p_main[[1]]/p_main[[2]]) 
p_out
#ggsave("viral_kinetic_compare.pdf",p_out,height=6,width=9)
#ggsave("viral_kinetic_hk.pdf",p_main[[1]],height=4,width=7)
##
#####

## end of script

#####