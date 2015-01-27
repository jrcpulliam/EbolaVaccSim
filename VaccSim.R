library(dplyr); library(data.table)
if(grepl('bellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
## Simulate SWCT vs RCT vs CRCT for SL
source('simFuns.R')

nGroups <- 16
nPerGroup <- 500
N <- nGroups*nPerGroup
hazPerMonth <- rgamma(nGroups, 1, 1)

idat <- data.table(id = 1:N, group = rep(1:nGroups, each = nPerGroup), vacc = 0, dis = 0, mort = 0)
idat

