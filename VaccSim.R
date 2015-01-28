library(dplyr); library(data.table)
if(grepl('bellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
## Simulate SWCT vs RCT vs CRCT for SL
source('simFuns.R')

yearToMonth <- 12/365.25
meanHaz <- .2 * yearToMonth  ## converting from per year
varC <- meanHaz^2/2
varI <- meanHaz^2/8

samp <- reParmRgamma(10^5, meanHaz, varC)
breaks <- 100
hist(samp, col = 'black', xlab = 'hazard / month', main = 'distribution of cluster hazards', 
     xlim = c(0, .05), breaks=breaks)
title(main=paste('average hazard = ', signif(meanHaz,2)), line = -2)

test <- makePop(20,300)
test <- setSW(test)
test <- setHazs(test, mean  = meanHaz, varClus = varC, varIndiv = varI)
test <- simSWCT(test, vaccEff = .6)

stest <- makeSurvDat(test)
t1 <- truncSurvDat(stest, 6)

summarise(group_by(t1, vacc), sum(infected))
## summarise(group_by(t1, vacc, cluster), sum(infected))

hist(test[,infectTime], breaks = seq(0, max(test[,infectTime])+1, by = 1), xlim = c(0, 4), xlab = 'months')
sum(test$infectTime < 30*4)


nGroups <- 16
nPerGroup <- 500
N <- nGroups*nPerGroup
hazPerMonth <- rgamma(nGroups, 1, 1)

idat <- data.table(id = 1:N, group = rep(1:nGroups, each = nPerGroup), vacc = 0, dis = 0, mort = 0)
idat

