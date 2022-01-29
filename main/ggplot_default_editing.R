require(ggplot2)

ggplot_theme <-     
        
        theme_set(
                
                theme_bw() +
                        
                        theme(axis.text.x = element_text(size = 4.6,
                                                         colour = "black"),
                              axis.text.y = element_text(size = 5,
                                                         colour = "black"),
                              panel.grid.major = element_blank(),
                              panel.background = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_blank(),
                              axis.line = element_line(colour = "black",size=.2),
                              axis.ticks = element_line(color = "black",size=.2),
                              axis.ticks.length = unit(.03,"cm"),
                              plot.margin = unit(c(0.1,0.2,0.01,0.1), "cm"),
                              axis.title.x = element_text(size = 5, face = 'bold',colour = "black",),
                              axis.title.y = element_text(size = 5, face = 'bold',colour = "black"),
                              plot.title = element_text(size = 7, face = 'bold',colour = "black"),
                              legend.text = element_text(size = 5, colour = "black"),
                              legend.title = element_text(size = 5, colour = "black"),
                              legend.key.size = unit(.3, 'cm'), #change legend key size
                              legend.key.height = unit(.25, 'cm'), #change legend key height
                              legend.key.width = unit(.3, 'cm'),
                              plot.title.position = 'plot')
                
        )
#
#


