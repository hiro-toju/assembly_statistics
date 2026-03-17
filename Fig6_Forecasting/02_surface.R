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
install.packages("devtools")
install.packages("vegan")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("cowplot")
install.packages("extrafont")
install.packages("RColorBrewer")
install.packages("scales")
install.packages("remotes")
install.packages("tidyverse")
install.packages("mgcv")
install.packages("plotly")
devtools::install_github("hiroakif93/R-functions/AnalysisHelper", force=TRUE)

library(AnalysisHelper)
## -- Loading Function and Library
load.lib( 'vegan', 'ggplot2', 'tidyverse', 'cowplot', 'extrafont','RColorBrewer', 'scales')

# -- Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)

# -- Load data table
sml <- readRDS('table/table_ela.rds')

############################################################################

i=2
	
	smlsub <- sml[sml$treat2=="OSN",]
	#smlsub <- smlsub[ order(smlsub$replicate.id), ]
	smlsub$color <-  rep(hcl.colors(110, "Spectral"), 8)
	## ================================================= ##
	## -- Surface slope estimating
	library(mgcv)
	library(plotly)
	
	## -- SS.Energy
	kk=8
	mod <- gam(Energy ~ te(rel.MDS1, k=5) + te(rel.MDS2, k=5) + ti(rel.MDS1, rel.MDS2, k=5), data=smlsub)
	
	mds1.seq <- seq(min(smlsub$rel.MDS1, na.rm=TRUE), max(smlsub$rel.MDS1, na.rm=TRUE), length=70)
	mds2.seq <- seq(min(smlsub$rel.MDS2, na.rm=TRUE), max(smlsub$rel.MDS2, na.rm=TRUE), length=70)
	
	predfun <- function(x,y){
	    newdat <- data.frame(rel.MDS1 = x, rel.MDS2=y)
	    predict(mod, newdata=newdat)
	}
	fit2 <- outer(mds1.seq, mds2.seq, Vectorize(predfun))
	
	cs <- scales::rescale(quantile(fit2, probs=seq(0,1,0.25)), to=c(0,1))
	frame2 <- plot_ly(data=smlsub, x = ~rel.MDS1, y= ~rel.MDS2, z= ~Energy) %>% 
	  add_trace(data=smlsub, x = ~rel.MDS1, y= ~rel.MDS2, z= ~Energy,
	            type = "scatter3d", mode = "markers",
	            marker = list(color = ~as.numeric(smlsub$time),
	                          colorscale = list(seq(0,1, length=8), 
	                                            c('firebrick', 'orange', 'khaki', 'greenyellow', 
	                                              'limegreen', 'dodgerblue', 'darkslateblue', 'darkslateblue')),
	                          #color ='orange', 
	                          size=5, legendgrouptitle=list(text='Energy', font='Arial'),
	                          line=list(width=1,color='black')),
	            opacity = 1)  %>% 
	  add_trace(data=smlsub[smlsub$time==1,], x = ~rel.MDS1, y= ~rel.MDS2, z= ~Energy,
	            text =c(1:8), mode = "text",
	            textfont = list(size=60, color='red'),
	            opacity = 1) %>%   
	  add_trace(data=smlsub[smlsub$time==110,], x = ~rel.MDS1, y= ~rel.MDS2, z= ~Energy,
	            text =c(1:8), mode = "text",
	            textfont = list(size=60, color='blue'),
	            opacity = 1) %>%           
	  add_trace(data=smlsub, x = ~mds1.seq, y= ~mds2.seq, z= ~t(fit2),
	            type = "surface", showscale = T,
	            hidesurface=F, opacity =0.7,
	            colorscale = list(cs, c('blue','lightblue','slategray', "tan", "indianred") ),
	            contours = list(z=list(show = TRUE, start = min(t(fit2)), end = max(t(fit2)), 
	                                   usecolormap=TRUE,  size=0.7, width=3)) )%>% 
	  layout( scene = list(xaxis = list(title = '', showticklabels=F,nticks=10, linewidth=7, gridwidth=3),
	                       yaxis = list(title = '', showticklabels=F, nticks=10, linewidth=7, gridwidth =3),
	                       zaxis = list(title = '', showticklabels=F, nticks=10, linewidth=7, gridwidth =3),
	                       aspectratio = list(x = .9, y = .9, z = 0.4),
	                       font='Arial') )	
	
	
	## ================================================= ##
	## -- Plotly
	
	mod2 <- gam(Energy ~ s(rel.MDS1, rel.MDS2), data= smlsub)
	
	rainbow(110)[smlsub$time]
	
	steps <- 40
	rel.MDS1 <- with(smlsub, seq(min(rel.MDS1), max(rel.MDS1), length = steps) )
	rel.MDS2 <- with(smlsub, seq(min(rel.MDS2), max(rel.MDS2), length = steps) )
	newdat <- expand.grid(rel.MDS1 = rel.MDS1, rel.MDS2 = rel.MDS2)
	Energy <- matrix(predict(mod2, newdat), steps, steps)
	
	perspcol = colorRampPalette(rev(brewer.pal(11,"RdBu")[-c(1:2, 4:5)]))
	# height of facets
	z.facet.center <- (Energy[-1, -1] + Energy[-1, -ncol(Energy)] + Energy[-nrow(Energy), -1] + Energy[-nrow(Energy), -ncol(Energy)])/4
	# Range of the facet center on a 100-scale (number of colors)
	range_z = 100
	z.facet.range<-cut(z.facet.center, range_z)

	pdf("result/soil_mediaA_point_col_time_v1.pdf", w=8, h=6, )
	p <- persp(rel.MDS1, rel.MDS2, Energy, theta = 30, phi=60, col = perspcol(range_z)[z.facet.range], ticktype="detailed", zlim=c(min(Energy), max(Energy)), zlab="Energy",
	           xlab="nMDS 1", ylab="nMDS 2")
	obs <- with(smlsub, trans3d(rel.MDS1, rel.MDS2, Energy, p))
	pred <- with(smlsub, trans3d(rel.MDS1, rel.MDS2, fitted(mod2), p))
	points(obs, col = alpha(rainbow(110)[smlsub$time], 0.5), pch = 16, cex=0.8)
	text(rel.MDS1[which(smlsub$time==1)], rel.MDS2[which(smlsub$time==1)], labels=1:8)
	
	dev.off()
	
	pdf("result//soil_mediaA_point_col_time_v2.pdf", w=8, h=6, )
	p <- persp(rel.MDS1, rel.MDS2, Energy, theta = 0, phi=60, col = perspcol(range_z)[z.facet.range], ticktype="detailed", zlim=c(min(Energy), max(Energy)), zlab="Energy",
	           xlab="nMDS 1", ylab="nMDS 2")
	obs <- with(smlsub, trans3d(rel.MDS1, rel.MDS2, Energy, p))
	pred <- with(smlsub, trans3d(rel.MDS1, rel.MDS2, fitted(mod2), p))
	points(obs, col = alpha(rainbow(110)[smlsub$time], 0.5), pch = 16, cex=0.8)
	text(rel.MDS1[which(smlsub$time==1)], rel.MDS2[which(smlsub$time==1)], labels=1:8)
	
	dev.off()

	pdf("result/soil_mediaA_point_col_replicate_v1.pdf", w=8, h=6, )
	p <- persp(rel.MDS1, rel.MDS2, Energy, theta = 30, phi=60, col = perspcol(range_z)[z.facet.range], ticktype="detailed", zlim=c(min(Energy), max(Energy)), zlab="Energy",
	           xlab="nMDS 1", ylab="nMDS 2")
	obs <- with(smlsub, trans3d(rel.MDS1, rel.MDS2, Energy, p))
	pred <- with(smlsub, trans3d(rel.MDS1, rel.MDS2, fitted(mod2), p))
	points(obs, col = alpha(c('firebrick', 'orange', 'khaki', 'greenyellow', 
	                          'limegreen', 'dodgerblue', 'darkslateblue', 'darkslateblue')[smlsub$replicate.id], 0.7), pch = 16, cex=0.8)
	text(rel.MDS1[which(smlsub$time==1)], rel.MDS2[which(smlsub$time==1)], labels=1:8)
	dev.off()
	
	
	pdf("result/soil_mediaA_point_col_replicate_v2.pdf", w=8, h=6, )
	p <- persp(rel.MDS1, rel.MDS2, Energy, theta = 0, phi=60, col = perspcol(range_z)[z.facet.range], ticktype="detailed", zlim=c(min(Energy), max(Energy)), zlab="Energy",
	           xlab="nMDS 1", ylab="nMDS 2")
	obs <- with(smlsub, trans3d(rel.MDS1, rel.MDS2, Energy, p))
	pred <- with(smlsub, trans3d(rel.MDS1, rel.MDS2, fitted(mod2), p))
	points(obs, col = alpha(c('firebrick', 'orange', 'khaki', 'greenyellow', 
	                          'limegreen', 'dodgerblue', 'darkslateblue', 'darkslateblue')[smlsub$replicate.id], 0.7), pch = 16, cex=0.8)
	text(rel.MDS1[which(smlsub$time==1)], rel.MDS2[which(smlsub$time==1)], labels=1:8)
	dev.off()

	pdf("result/soil_mediaA_point_col_replicate_v2.pdf", w=8, h=6, )
	p <- persp(rel.MDS1, rel.MDS2, Energy, theta = 20, phi=40, col = perspcol(range_z)[z.facet.range], ticktype="detailed", zlim=c(min(Energy), max(Energy)), zlab="Energy",
	           xlab="nMDS 1", ylab="nMDS 2")
	obs <- with(smlsub, trans3d(rel.MDS1, rel.MDS2, Energy, p))
	pred <- with(smlsub, trans3d(rel.MDS1, rel.MDS2, fitted(mod2), p))
	points(obs, col = alpha(c('firebrick', 'orange', 'khaki', 'greenyellow', 
	                          'limegreen', 'dodgerblue', 'darkslateblue', 'darkslateblue')[smlsub$replicate.id], 0.7), pch = 16, cex=0.8)
	text(rel.MDS1[which(smlsub$time==1)], rel.MDS2[which(smlsub$time==1)], labels=1:8)
	dev.off()
	