if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
require(RColorBrewer); require(boot)
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R', 'extractFXN.R'), source)

thing <- 'Equip-rand'
out <- extractSims(thing, verb=0)
load(file=file.path('BigResults',paste0(thing, '.Rdata')))
finTrials[order(gs), list(tcalMean = mean(tcal), power = mean(vaccGood)), list(vaccEff, trial, gs, ord, delayUnit)]
finTrials[tcal<168, list(nsim=length(tcal), tcalMean = mean(tcal), power = mean(vaccGood)), list(vaccEff, trial, gs, ord, delayUnit)]
finInfo[cat=='allFinalEV', list(caseTot, simNum)]

pf <- makePowTable(finTrials, ver=0)
pf[order(gs), list(gs, trial, ord, delayUnit, vaccEff, biasNAR, cvrNAR, vaccGood)]
names(pf)

####################################################################################################
## extract counterfactuals
load('data/vaccProp1.Rdata')
vaccProp <- vaccProp1
vaccProp[, simNum:=1:length(vaccEff)]

thing <- 'Equip-randCFs'
fincfs <- extractCFs(thing, verb=0)
load(file=file.path('BigResults',paste0(thing, '.Rdata')))
fincfs <- fincfs[nbatch<161] ## next 160-320 are redundant
fincfs <- merge(fincfs, vaccProp, by = 'simNum', all.y=F) ## copy vaccProp in there (should do this in analysisFuns.R later

## ## Look at counterfactuals
## fincfs[, list(caseTot= mean(caseTot), caseTotmin= min(caseTot), caseTotmax= max(caseTot), n=length(caseTot)) , list(simNum, cf)]
## ## 
## fincfs[cf=='VRpop', list(caseTot= mean(caseTot), caseTotmin= min(caseTot), caseTotmax= max(caseTot), n=length(caseTot), vaccEff= unique(vaccEff)) , list(simNum)]
## summary(lm(caseTot ~ vaccEff, data = fincfs[cf=='VRpop']))
## ## 
## jpeg('Figures/test.jpg')
## p <- ggplot(fincfs[cf=='VRpop'], aes(vaccEff, caseTot, group=simNum, color=simNum))
## ## p + geom_point(aes(color=simNum))
## p + geom_boxplot()
## graphics.off()
## ## 
## jpeg('Figures/test2.jpg')
## hist(fincfs[cc==1,vaccEff])
## graphics.off()
## ## 
## fincfs[cf=='NTpop', var(caseTot)] ## far more variance between hazard simulations than within them
## fincfs[cf=='NTpop', var(caseTot), list(simNum)]

####################################################################################################
## compare to factuals

fincfs[, length(vaccEff), list(simNum,cf)]
finTrials[, length(vaccEff), simNum] ## 4 types of simulations
finTrials[, list(numSim = length(vaccEff)), list(simNum, gs, ord)] ## these are the groupings
finInfo[cat=='allFinalEV', length(caseTot), list(simNum, gs, ord)] ## 4 types of simulations
setkey(finInfo, gs, ord, simNum, cat)

class(fincfs$cc) <- 'numeric'
setkey(fincfs,  simNum, cc)
names(fincfs)

## Merge (big merger but worth it in speed later)
fin <- merge(finInfo, fincfs, by ='simNum', all.y=F, allow.cartesian = T)
setkey(fin, gs, ord.x, simNum, cf, cat)

## Get difference between caseTot
ctot <- fin[cat=='allFinalEV', list(caseTotvsCF = caseTot.x - caseTot.y), list(gs, ord.x, simNum, cf)]
ctot[, mean(caseTotvsCF), list(gs, ord.x, simNum, cf)]
