if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
source('simFuns.R')

parms <- makeParms()
samp <- reParmRgamma(10^5, parms$mu, parms$varClus)
breaks <- 100
hist(samp/monthToDays, col = 'black', xlab = 'hazard / month', main = 'distribution of cluster hazards', 
     xlim = c(0, .05), breaks=breaks)
title(main=paste('average hazard = ', signif(parms$mu/monthToDays,2)), line = -2)

swctSims <- simNtrials(parms=makeParms('SWCT'), N = 50)
rctSims <- simNtrials(parms=makeParms('RCT'), N = 50)

resSWCT <- simTrial(makeParms('SWCT')) #, numClus = 2, clusSize = 6))
resSWCT
resRCT <- simTrial(makeParms('RCT')) #, numClus = 2, clusSize = 6))
resRCT

mu <- makeParms()$mu

censSurvDat(resSWCT, 40)[,active]
censSurvDat(resRCT, 40)[,active]

resSWCT <- simTrial(makeParms('SWCT', varClus = mu^2/2, varIndiv = mu^2/1, numClus = 20, clusSize = 300, vaccEff = .6))
tt <- 60
for(tt in 7*(6:20)) print(doCoxPH(censSurvDat(resSWCT, tt)))
doGlmer(censSurvDat(resSWCT, tt),T)
doGlmer(censSurvDat(resSWCT, tt),F)

summTrial(censSurvDat(resRCT$st, 30))
doCoxPH(censSurvDat(resRCT$st, 70))


firstStop(simTrial(makeParms('SWCT')), verb=1)
firstStop(simTrial(makeParms('RCT')), verb=1)

t1 <- censSurvDat(resSWCT$st, 180)

doCoxPH(censSurvDat(stest$st, 62/30))
summarise(group_by(t1, vacc), sum(infected))
## summarise(group_by(t1, vacc, cluster), sum(infected))

hist(test[,infectTimeTrunc], breaks = 0:12, xlim = c(0, 12), col = 'black', xlab = 'infection times (months)')



nGroups <- 16
nPerGroup <- 500
N <- nGroups*nPerGroup
hazPerMonth <- rgamma(nGroups, 1, 1)

idat <- data.table(id = 1:N, group = rep(1:nGroups, each = nPerGroup), vacc = 0, dis = 0, mort = 0)
idat

