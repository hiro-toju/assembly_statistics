set.seed(123)
n.core=12

library(ggplot2); library(vegan);library(gridExtra);library(lemon);library(tidyr);library(markdown)
library(stringr);library(ggforce);library(cowplot);library(ggnewscale); library(reshape2);library(ggtext)
library(ecoregime)
source("./functions.R")

library("doParallel")
library('tidyverse')
library('plyr')
library('gtools')
library('ggsci')
library('igraph')
library('tidygraph')
library('RColorBrewer')
library("stringdist")
library('vegan')
library("plotly")
library("parallel")
library("dplyr")
library("mgcv")
library("rELA")
library(parallel)
library(foreach)
library(doParallel)

dir.create("ELA_OUT")

# -- read data from previous dir
exp <- read.csv("table/expdata.csv")
data <- readRDS("table/List_data_rep.rds")

tax <- c("cla","ord","fam","gen","ASV")
Medium <- c("G","GL","GC","GLC","HG","HGL","HGC","HGLC")

#### -- Data -- ####
#- ndata[[i]][[j]]:
#- i = the number of taxonomy levels, j = the number of Medium

GGG <- function(data, M=Medium){ #- data <- data[[2]]
  data <- data[,!(colnames(data) %in% c("Sample_name", "Plate"))]
  data$Day <- (data$Day)/2
  
  ans <- list()
  for(i in 1:length(M)){
    ans[[i]] <- data[(data$Medium == M[i]),]
    ans[[i]]$Location <- as.integer(factor(ans[[i]]$Location, levels=unique(ans[[i]]$Location)))
  }
  
  return(ans)
}

ndata <- list()
for(i in 1:length(tax)){
  ndata[[i]] <- GGG(data=data[[i]])
}

saveRDS(ndata, file="table/Hayashi_ndata.rds")

####################################################################
####################################################################
# selecting taxonomic level and medium

x <- "asvGLC"
input <- ndata[[5]][[4]]

####################################################################
####################################################################
###--- RETRA-EDR function
Msegs <- 6 #-- Calculation in multiple Minsegs

trajectories <-  input$Location
states <-  as.integer(input$Day)

dist <- vegdist(input[,-c(1:3)], method="bray")

edrout <- retra_edr(d=dist, trajectories=trajectories, states=states, minSegs=6)

srn <- 1; srt <- "T1"
for(t in 2:length(edrout)){
  if(edrout[[srn]]$Length < edrout[[t]]$Length){
    srn <- t; srt <- sprintf("T%s",as.character(t))
  }
}

####################################################################
####################################################################
# Comparison between cmdscale and EDR

pdf(sprintf("ELA_OUT/cmdscale_vs_EDR_%s.pdf", x), height=5, width=9)
par(mfrow = c(1, 2))

plot(cmdscale(dist), col=input$Day, pch=input$Day)

plot(x = edrout, d = dist, trajectories = trajectories, states = states,
     select_RT = srt, #if we need to highlight some RT
     traj.colors = "lightblue", RT.colors = "orange", sel.color = "darkgreen",
     link.lty = 1, asp = 1)

dev.off()
####################################################################
####################################################################

#########################################################################
## ELA

relmat <- input[,-c(1:3)]
envmat <- input[,c(1:2)]

rel <- 0.005
occ <- 0.05
prun <- 0.1
qth <- 10^-5

eladata <- Formatting(relmat, envmat, normalize=1, c(rel, occ, 1 - occ), grouping=0, grouping_th=NULL)
ocmat <- eladata[[1]]
enmat <- eladata[[3]]

saveRDS(ocmat, file=sprintf("ELA_OUT/ocmat_asvG_%s_%s.rds", rel, occ))
write.csv(x = ocmat, file=sprintf("ELA_OUT/ocmat_asvG_%s_%s.csv",  rel, occ)) 

dim(ocmat)
sum(ocmat)

sa <- runSA(ocmat=as.matrix(ocmat), enmat=NULL, qth=qth, rep=256, threads=n.core)
saveRDS(sa, file=sprintf("ELA_OUT/runSA_%s_%s_%s.rds", rel, occ, x))
#sa <- readRDS(file=sprintf("ELA_OUT/runSA_%s_%s_%s.rds", rel, occ, x))

gstb <- gStability(sa, ocmat, enmat=NULL, th=prun, threads=n.core)
saveRDS(gstb, file=sprintf("ELA_OUT/gStability_%s_%s_%s.rds", rel, occ, x))
#gstb <- readRDS(file=sprintf("ELA_OUT/gStability_%s_%s_%s.rds", rel, occ, x))


############################################################################
# Bray-Curtis PCoA
############################################################################
############################################################################
# Preparing data
pcoa <- cmdscale(dist)
colnames(pcoa) <- c("PCoA.1", "PCoA.2")

############################################################################
# Preparing data

df <- cbind(pcoa, gstb[[1]])
sslist <- gstb[[4]][[4]][[2]]
colnames(sslist)[1] <- "stable.state.id"
#df <- left_join(df, sslist, by="stable.state.id")

############################################################################
# Graphic parameters

n.col <- length(unique(df$stable.state.id))
colvec <- c(brewer.pal(n = n.col, name = "Set1"))
col1 <- colvec[1:n.col]

par(oma = c(0.5, 0.5, 0.5, 0.5)) 
par(mar = c(0.5, 0.5, 0.5, 0.5)) 

#######################################################
# Landscape plotting

mod1 <- gam(e.realize ~ s(PCoA.1, PCoA.2), data= df)

##
mds1.seq <- seq(min(df$PCoA.1, na.rm=TRUE), max(df$PCoA.1, na.rm=TRUE), length=70)
mds2.seq <- seq(min(df$PCoA.2, na.rm=TRUE), max(df$PCoA.2, na.rm=TRUE), length=70)

predfun <- function(x,y){
  newdat <- data.frame(PCoA.1 = x, PCoA.2=y)
  predict(mod1, newdata=newdat)
}
fit <- outer(mds1.seq, mds2.seq, Vectorize(predfun))

###
## -- Plotly
cs <- scales::rescale(quantile(fit, probs=seq(0,1,0.25)), to=c(0,1))

names(cs) <-NULL
#df$color=colors
frame6 <- plot_ly(data=df, x = ~PCoA.1, y= ~PCoA.2, z= ~e.realize) %>% 
  add_trace(data=df, x = ~PCoA.1, y= ~PCoA.2, z= ~e.realize,
            type = "scatter3d", mode = "markers",
            marker = list(color = ~(envmat$Day),
                          colorscale = list(seq(0,1, length=6), 
                                            c('firebrick', 'orange', 'khaki', 
                                              'limegreen', 'dodgerblue', 'darkslateblue')),
                          size=2, legendgrouptitle=list(text='e.realize', font='Arial'),
                          line=list(width=1,color='black')), opacity = 1)  %>% 
  add_trace(data=df, x = ~mds1.seq, y= ~mds2.seq, z= ~t(fit),
            type = "surface", showscale = F,
            hidesurface=F, opacity =0.7,
            colorscale = list(cs, c('blue','lightblue','slategray', "tan", "indianred") ),
            contours = list(z=list(show = TRUE, start = min(t(fit)), end = max(t(fit)), 
                                   usecolormap=TRUE,  size=0.7, width=3))) %>% 
  layout( scene = list(xaxis = list(title = '', showticklabels=F,nticks=10, linewidth=7, gridwidth=3),
                       yaxis = list(title = '', showticklabels=F, nticks=10, linewidth=7, gridwidth =3),
                       zaxis = list(title = '', showticklabels=F, nticks=10, linewidth=7, gridwidth =3),
                       aspectratio = list(x = .9, y = .9, z = 0.4),
                       font='Arial') )	

frame6


####################################################################



