library(dplyr); library(data.table)
if(grepl('bellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
## Simulate SWCT vs RCT vs CRCT for SL
source('simFuns.R')

yearToMonth <- 12/365.25
meanHaz <- .2 * yearToMonth  ## converting from per year
varC <- .1 * yearToMonth^2   ## cluster level
varI <- .05 * yearToMonth^2  ## individual level

breaks <- seq(0,.2, by = .01)
hist(reParmRgamma(10^5, meanHaz, varC), col = 'black', xlab = 'hazard / month', main = 'distribution of cluster hazards', 
     xlim = c(0, .1), breaks=breaks)
title(main=paste('average hazard = ', signif(meanHaz,2)), line = -2)

test <- makePop(20,300)
test <- setSW(test)
test <- setHazs(test, mean  = meanHaz, varClus = varC, varIndiv = varI)
test <- simSWCT(test, vaccEff = .6)

stest <- makeSurvDat(test)

t1 <- truncSurvDat(stest, 6)

summarise(group_by(t1, vacc), sum(infected))
## summarise(group_by(t1, vacc, cluster), sum(infected))

mod <- coxph(Surv(startTime, endTime, infected) ~ vacc, data=t1)
modF <- coxph(Surv(startTime, endTime, infected) ~ 
             vacc + frailty.gamma(cluster, eps=1e-10, method="em", sparse=0),
             outer.max=1000, iter.max=10000,
             data=t1)
summary(mod)
summary(modF)

hist(test[,infectTime], breaks = seq(0, max(test[,infectTime])+1, by = 1), xlim = c(0, 4), xlab = 'months')
sum(test$infectTime < 30*4)


nGroups <- 16
nPerGroup <- 500
N <- nGroups*nPerGroup
hazPerMonth <- rgamma(nGroups, 1, 1)

idat <- data.table(id = 1:N, group = rep(1:nGroups, each = nPerGroup), vacc = 0, dis = 0, mort = 0)
idat

