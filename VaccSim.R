if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R'), source)

parms <- makeParms()
samp <- reParmRgamma(10^5, parms$mu, parms$varClus)
breaks <- 100
hist(samp/monthToDays, col = 'black', xlab = 'hazard / month', main = 'distribution of cluster hazards', 
     xlim = c(0, .05), breaks=breaks)
title(main=paste('average hazard = ', signif(parms$mu/monthToDays,2)), line = -2)

swctSims <- simNtrials(parms=makeParms('SWCT'), N = 3)
rctSims <- simNtrials(parms=makeParms('RCT'), N = 3)

resSWCT <- simTrial(makeParms('SWCT')) #, numClus = 2, clusSize = 6))
resSWCT
resRCT <- simTrial(makeParms('RCT')) #, numClus = 2, clusSize = 6))
head(resRCT$pop,50)

mu <- makeParms()$mu

censSurvDat(resSWCT, 40)[,active]
censSurvDat(resRCT, 40)[,active]

resSWCT <- simTrial(makeParms('SWCT', numClus = 20, clusSize = 500, vaccEff = .6)) ## varClus = mu^2/2, varIndiv = mu^2/1, 
tt <- 200
temp <- censSurvDat(resSWCT, tt)
summTrial(temp)
#for(tt in 7*(6:20)) 
doCoxPH(temp,'coxme')
doGlmer(temp,T)
doGlmer(temp,F)

## not an issue with 0s it seems
list(summarise(group_by(temp, cluster), sum(infected))
                               , summarise(group_by(temp, cluster, vacc), sum(infected))
                               , summarise(group_by(temp, vacc), sum(infected))
                               )

summTrial(censSurvDat(resRCT$st, 30))
doCoxPH(censSurvDat(resRCT$st, 70))

firstStop(simTrial(makeParms('SWCT')), verb=0)
seqStop(simTrial(makeParms('SWCT')), verb=0)

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

