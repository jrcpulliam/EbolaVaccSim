if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
require(RColorBrewer); require(boot); require(data.table); require(vioplot)
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
## 
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
setkey(finTrials, gs, ord, simNum)

finTrials <- finTrials[trial!='SWCT']
finInfo <- finInfo[trial!='SWCT']

## merge them so we have power & speed info here too
finit <- finInfo[finTrials[, list(gs, ord, simNum, tcal, vaccCases, contCases, vaccGood, vaccBad, vaccEff, PHU)]]


class(fincfs$cc) <- 'numeric'
setkey(fincfs,  simNum, cc)
names(fincfs) 

finit[cat=='allFinalEV'][1:2]

breaks <- seq(0,200, by = 5)
pdf('Figures/tp1.pdf')
par(mfrow=c(2,1))
finit[cat=='allFinalEV' & gs==T, 
hist(tcal, xlab = 'DALYS lost', freq = F, breaks = breaks, col = 'black', bty = 'n', xlim = range(breaks), plot=T, main=paste(ord)), 
list(ord)]
dev.off()

finit[cat=='allFinalEV' & gs==T, mean(tcal<168), ord]
finit[cat=='allFinalEV', mean(caseTot), list(ord, gs)]

pdf('Figures/tp1.pdf')
par(mfrow=c(2,1))
for(cc in c('allFinalEV', 'allFinal_noEV')) {
    finit[, lab:=factor(paste0(c('','gs-')[as.numeric(gs==T)+1],ord))]
    finit[, col:=rainbow(4)[as.numeric(lab)]]
    plot(0,0, type = 'n', xlab = 'DALYS lost', col = 'black', bty = 'n', xlim = c(0,200), ylab='',ylim = c(0, 0.02), las = 1, main =cc)
    finit[cat==cc, lines(density(caseTot), col = col[1] , lwd = 2), list(ord, gs)]
}
finit[, legend('topright', leg = unique(lab), col = unique(col), lwd = 2, bty = 'n')]
dev.off()



## Merge (big merger but worth it in speed later)
fin <- merge(finIT, fincfs, by ='simNum', all.y=F, allow.cartesian = T)
class(fin$simNum) <- class(fin$cc)<- 'integer'
setkey(fin, gs, ord.x, simNum, cf, cat)

## Get difference between caseTot
ctot <- fin[cat %in% c('allFinalEV', 'allFinal_noEV'), 
            list(cc, vaccGood, vaccBad, tcal, caseTotvsCF = caseTot.x - caseTot.y, caseTotF = caseTot.x, caseTotCF = caseTot.y), 
            list(gs, ord.x, simNum, cf, cat)]
setkey(ctot, simNum, cc, cf)
satr <- CJ(simNum = 1:2080, cc = 1:500, cf = c('NTpop','VRpop'))

## Merge to make sure that empty columns exist
ctot <- ctot[satr, allow.cartesian=T]
ctot


## Get difference between caseTot
ctot <- fin[cat %in% c('allFinalEV', 'allFinal_noEV'),  ## [1] below is because rows are replicated by cf, except for cf & caseTot.y & SAE
            list(cc[1], vaccGood[1], vaccBad[1], tcal[1],  caseTotF = caseTot.x[1], caseTotNT = caseTot.y[cf=='NTpop'], caseTotVR = caseTot.y[cf=='VRpop']), 
            list(gs, ord.x, simNum, cat)]
ctot

ctot[, list(ctF = mean(caseTotF), ctCF = mean(caseTotCF), diff = mean(caseTotvsCF)), list(gs, ord.x, simNum, cf, cat)]

cs <- ctot[order(cat), list(ctF = mean(caseTotF), ctCF = mean(caseTotCF), diff = mean(caseTotvsCF), 
                            power = mean(vaccGood), tcal = mean(tcal)), 
           list(gs, ord.x, cf, cat)]
save(cs, file = file.path('Results','cs.Rdata'))

load(file = file.path('Results','cs.Rdata'))
class(cs$gs) <- 'logical'
setkey(cs, gs, ord.x)

par(mfrow=c(2,2))
cs[, {


}
 , list(gs, ord.x)] 

ctot

cs[, lab:=paste0(c('gs','')[as.numeric(gs) + 1], 'RCT-',c('none','TU')[as.numeric(ord.x=='TU') + 1])]
cs[, lab:=factor(lab)]

pdf('Figures/test.pdf')
plot(0,0, type = 'n', bty = 'n', ylim = c(0,.5), xlim = cs[.(T, 'TU'), range(ctF,ctCF)], las = 1, xlab='cases', ylab = 'power')
cs[cf=='NTpop' & cat=='allFinalEV', points(ctF, power, pch = 15, cex = 2, col = as.numeric(lab))]
legend('bottom', leg = cs[, unique(lab)], col = cs[, as.numeric(unique(lab))], bty = 'n', pch = 15)
abline(v = cs[cf=='NTpop', ctCF], lty=2)
abline(v = cs[cf=='VRpop', ctCF], lty=2)
dev.off()


pdf('Figures/test.pdf')
plot(0,0, type = 'n', bty = 'n', ylim = c(0,.5), xlim = cs[.(T, 'TU'), range(ctF,ctCF)], las = 1, xlab='cases', ylab = 'power')
cs[cf=='NTpop' & cat=='allFinalEV', points(ctF, power, pch = 15, cex = 2, col = as.numeric(lab))]
legend('bottom', leg = cs[, unique(lab)], col = cs[, as.numeric(unique(lab))], bty = 'n', pch = 15)
abline(v = cs[cf=='NTpop', ctCF], lty=2)
abline(v = cs[cf=='VRpop', ctCF], lty=2)
dev.off()

cs[cf=='NTpop' & cat=='allFinalEV']

pdf('Figures/test.pdf')
ctot

