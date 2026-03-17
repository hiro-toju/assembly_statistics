############################################################################
####
#### R script for Fujita (2022)
####
#### Removing outlier from data
#### 2022.03.16 Fujita
#### 2022.06.23 Toju
#### 2023.12.21 Toju
#### R 4.3.2
#### 
############################################################################
set.seed(123)

library("Rcpp")
library("RcppArmadillo")
library("doParallel")
library('tidyverse')
library('gsubfn')
library('zoo')
library('plyr')
library('gtools')
library('ggsci')
library('igraph')
library('tidygraph')
library('RColorBrewer')
library("stringdist")
library("purrr")
library("plot3D")
library("rELA")
library("foreach")
library("parallel")
library("dplyr")
#source("showGELA3D_mod_240808.R")
source("gelsobj2.R")

n.core=12
#########################################################################

# loading data
baseabtable <- readRDS('table/ELA_matrix_each_fun_Family.rds')
env <- readRDS('table/env_fun_1664.rds')
basemetadata <- data.frame(apply(env, 2, scale))
rownames(basemetadata) <- rownames(env)

#########################################################################
#########################################################################
## Analysis with environmental data (metadata)
## Setting 1

rel <- 0.001
occ <- 0.10
prun <- 0.1
qth <- 10^-5

#list[ocmat, abmat, enmat, samplelabel, specieslabel, factorlabel] <- Formatting(baseabtable, basemetadata, 1, c(rel, occ, 1 - occ))
input <- Formatting(baseabtable, basemetadata, normalize=1, c(rel, occ, 1 - occ), grouping=0, grouping_th=NULL)
ocmat <- input[[1]]
enmat <- input[[3]]

saveRDS(ocmat, file=sprintf("result/ocmat_fun_Family_%s_%s.rds", rel, occ))
write.csv(x = ocmat, file=sprintf("result/ocmat_fun_Family_%s_%s.csv",  rel, occ)) 

dim(ocmat)
sum(ocmat)

#sa <- runSA(ocmat=as.matrix(ocmat), enmat=as.matrix(enmat), rep=256, threads=n.core)
#sa <- runSA(ocmat=as.matrix(ocmat), enmat=as.matrix(enmat), qth=qth, rep=256, threads=n.core)
#saveRDS(sa, file=sprintf("table/runSA_fun_Family_%s_%s.rds", rel, occ))

sa <- readRDS(file=sprintf("table/runSA_fun_Family_%s_%s.rds", rel, occ))

#########################################################################
## GStability Manual

#gStability returns a list of 4 elements: the first two are the dataframe for pruned/non-pruned energy landscape, respectively. In addition to the dataframe of Stability it includes e.tipping (energy of tipping point) and energy.barrier (height of energy from observed state to the tipping point).The third output is a list of parameters (h, g, j, h+g*env) and a summary table of stable states, and the fourth output is a list encapsulating the inputs required for the various plots.

#output of gStability:
#[[1]]: data.frame(energy.gap, ss.entropy, energy.barrier, e.realize, e.stable, e.tipping, state.id, stable.state.id)
#[[2]]: data.frame(energy.gap.np, ss.entropy.np, energy.barrier.np, e.realize, e.stable.np, e.tipping.np, state.id.np, stable.state.id.np)
#[[3]]: w/ enmat: list(list(list(he, je, ge, hge), data.frame(sstable)), ...); w/o enmat: list(list(he, je, ge, hge), data.frame(sstable))
#[[4]]: w/ enmat: list(list(ocmat, env, sa, ela, elanp), ...); w/o enmat: list(ocmat, env, sa, ela, elanp)

#########################################################################
# gStability

# assuming mean values of environmental variables
#gstb <- gStability(sa, ocmat, enmat=NULL, th=prun, threads=n.core) # assuming mean values of environmental variables
#saveRDS(gstb, file=sprintf("table/gStability_meanEnv_fun_Family_%s_%s.rds", rel, occ))
#write.csv(x = gstb[[1]], file=sprintf("result/gStability_meanEnv_fun_Family_%s_%s_prun%s.csv",  rel, occ, prun)) #pruned
#write.csv(x = gstb[[2]], file=sprintf("result/gStability_meanEnv_fun_Family_%s_%s.csv",  rel, occ)) #non-pruned

gstab <- readRDS(file=sprintf("table/gStability_meanEnv_fun_Family_%s_%s.rds", rel, occ))

#########################################################################
# PCA 

# if enmat=NULL, remove "[[sample.id]]"
ocmat <- gstb[[4]][[1]]
env <- gstb[[4]][[2]]
sa <- gstb[[4]][[3]]
ela <- gstb[[4]][[4]]

pdf(file=sprintf("result/PCplot_fun_Family_%s_%s.pdf", rel, occ))
PCplot(ocmat, sa, ssrep=ela[[2]], pruned=FALSE)
dev.off()

pdf(file=sprintf("result/PCplot_fun_Family_%s_%s_prun%s.pdf", rel, occ, prun))
PCplot(ocmat, sa, ssrep=ela[[2]])
dev.off()

#########################################################################

elanp <- gstb[[4]][[5]]

pdf(file=sprintf("result/showDG_gstab5_env_fun_Family_%s_%s.pdf", rel, occ))
showDG(elanp, ocmat, "not-pruned")
dev.off()

pdf(file=sprintf("result/showDG_gstab4_env_fun_Family_%s_%s_prun%s.pdf", rel, occ, prun))
showDG(ela[[1]], ocmat, "pruned")
dev.off()

#########################################################################
#########################################################################
# ELA gradient
#########################################################################

# Available P

gela <- GradELA(sa=sa, eid="available_P", # Specify the label or position of an environmental factor
                enmat=enmat, env=NULL, range=NULL, steps=32, th=prun, threads=n.core) #[[1]]: return value of ELA function for each step, [[2]]: value of environmental factor for each step, [[3]]: specified environmental factor

saveRDS(gela, file=sprintf("table/gela_meanEnv_fun_Family_availableP_%s_%s_prun%s.rds", rel, occ, prun))

gela <- readRDS(file=sprintf("table/gela_meanEnv_fun_Family_availableP_%s_%s_prun%s.rds", rel, occ, prun))


gelsobj <- GELSObj2(gela, sa, threads=n.core)
pdf(file=sprintf("result/ELAgradient_fun_Family_availableP_%s_%s_prun%s_mod.pdf",  rel, occ, prun))
showGELA3D(gelsobj)
dev.off()

#########################################################################

