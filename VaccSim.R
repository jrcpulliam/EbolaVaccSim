if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

p1 <- simTrial(makeParms('CRCT',small=F, ord='TU'))
s1 <- makeSurvDat(p1,'pop',T)
s1 <- activeFXN(s1,'st',T)
sc1 <- censSurvDat(s1, 49)
summTrial(censSurvDat(s1, 365,'st'))
summTrial(censSurvDat(s1, 365,'stActive'))
doCoxPH(censSurvDat(s1, 154),br=F)
t1 <- seqStop(s1, verbose=0)
t1 <- endT(t1)
t1 <- makeCaseSummary(t1)

summTrial(t1$st)

## pairs matched for randomization (if matching)
browser()
popH$pair <- popH[, cluster %% (numClus/2)]
popH[pair==0, pair:=numClus/2] 

swctSims <- simNtrials(parms=makeParms('SWCT'), N = 10)
rctSims <- simNtrials(parms=makeParms('RCT'), N = 10)
crctSims <- simNtrials(parms=makeParms('CRCT'), N = 10)

resSWCT <- simTrial(makeParms('SWCT')) #, numClus = 2, clusSize = 6))
resSWCT
resRCT <- simTrial(makeParms('RCT')) #, numClus = 2, clusSize = 6))
resCRCT <- simTrial(makeParms('CRCT')) #, numClus = 2, clusSize = 6))
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

firstStop(simTrial(makeParms('SWCT')), verb=1)
seqStop(simTrial(makeParms('SWCT')))
seqStop(simTrial(makeParms('RCT')))
seqStop(simTrial(makeParms('CRCT')))

firstStop(simTrial(makeParms('RCT')), verb=0)

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

## Temporal variation in incidence every week *
## reording of vaccination sequence in time *
## Vaccination upon trial finishing in RCT/CRCT*
## calculate # infected by active & by maxInfectDays*
## Test false positives
## add CFR
## Equipoise calculations

## CRCT matched means you have less control groups active early on
