if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
require(RColorBrewer); require(boot)
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R', 'extractFXN.R'), source)

thing <- 'Equip1'
finTrials <- extractSims(thing, verb=0)
finTrials[order(gs), list(tcalMean = mean(tcal), power = mean(vaccGood)), list(vaccEff, trial, gs, ord, delayUnit)]
finTrials[tcal<168, list(nsim=length(tcal), tcalMean = mean(tcal), power = mean(vaccGood)), list(vaccEff, trial, gs, ord, delayUnit)]
load(file=file.path('BigResults',paste0(thing, '.Rdata')))

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

## Look at counterfactuals
fincfs[, list(caseTot= mean(caseTot), caseTotmin= min(caseTot), caseTotmax= max(caseTot), n=length(caseTot)) , list(simNum, cf)]

fincfs[cf=='VRpop', list(caseTot= mean(caseTot), caseTotmin= min(caseTot), caseTotmax= max(caseTot), n=length(caseTot), vaccEff= unique(vaccEff)) , list(simNum)]

summary(lm(caseTot ~ vaccEff, data = fincfs[cf=='VRpop']))

jpeg('Figures/test.jpg')
p <- ggplot(fincfs[cf=='VRpop'], aes(vaccEff, caseTot, group=simNum, color=simNum))
## p + geom_point(aes(color=simNum))
p + geom_boxplot()
graphics.off()

fincfs[cf=='NTpop', var(caseTot)] ## far more variance between hazard simulations than within them
fincfs[cf=='NTpop', var(caseTot), list(simNum)]

summary(aov(caseTot ~ simNum, data = fincfs))

## compare to factuals
