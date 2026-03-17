############################################################################
####
#### R script for Fujita (2019)
####
#### Visualizing community assembly
#### 2019.12.09 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS_6/')
#### 
############################################################################

## -- Loading Function and Library
install.packages("devtools")
install.packages("vegan")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("cowplot")
install.packages("extrafont")
install.packages("RColorBrewer")
install.packages("scales")
install.packages("remotes")
devtools::install_github("hiroakif93/R-functions/AnalysisHelper", force=TRUE)

library(AnalysisHelper)
load.lib( 'vegan', 'ggplot2', 'tidyr', 'cowplot', 'extrafont','RColorBrewer', 'scales')

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load data table
dlist  <- readRDS("table/table_assembly.rds")
sml <- readRDS('table/table_ela.rds')

taxa <- as.matrix(readRDS("table/table_taxa.rds"))
#color <- readRDS('table/99_color_palette_v2.rds')

############################################################################
glist1 <- c()
#for(i in names(dlist)[-7]){ #
l='Genus'
i=names(dlist)[2]
    
    ## ============================================================ ##
    ## -- Extracting one treatment matrix
    
    smlsub <- sml[sml$treat1==i, ]
    ts <- dlist[[i]]
    rel.ts <- ts/rowSums(ts)
    a <- smlsub[which(smlsub[,'abruptness_rel_tw5_tp1']>0.5), ]
    ## ======================================= ##
    ## -- Visualizing community assembly
    
    agg.ts <- Taxa.mat(ts[rownames(smlsub),], taxa, l)
    
    barobj = barplot_obj(agg.ts, specificName = "unidentified")
    
    lf <- gather(cbind(smlsub, barobj[[1]]), 
                 key, value, -c(1:(ncol(smlsub))))|>
      mutate(replicate.id = paste("Replicate", replicate.id))
    
    lf$key <- factor(lf$key, levels=rev(names(barobj[[2]]))) 
    
    g2 <- ggplot(lf)+
      geom_area(aes(x=as.numeric(time), y=value, fill=key), position='fill', color='grey30', stat='identity', size=0.1)+
      facet_wrap(~replicate.id,ncol=2, strip.position='top')+
      scale_fill_manual( values=barobj[[2]] ) + 
      labs(x= "Day", subtitle=i, y= "Relative abundance")+
      guides(fill=guide_legend(title='',ncol=1, reverse=FALSE))+
      theme_text(family="Arial", expandCx = TRUE, expandCy = TRUE, unitsize = 0.5)+
      theme(strip.text.x = element_text(hjust = 0)) 
    ggsave("result/community_assembly.pdf", g2, w=15,h=12, units = "cm", device = cairo_pdf)
