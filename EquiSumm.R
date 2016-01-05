if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T))
require(RColorBrewer); require(boot); require(data.table); require(vioplot)
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R', 'extractFXN.R'), source)

####################################################################################################
## extract factuals
thing <- 'Equip-rand3pit'
## out <- extractSims(thing, verb=0)
load(file=file.path('BigResults',paste0(thing, '.Rdata')))
finTrials[order(gs), list(tcalMean = mean(tcal), power = mean(vaccGood), length(tcal)),
          list(trial, gs, ord, delayUnit, propInTrial)]
finTrials[tcal<168 & vaccEff>0, list(nsim=length(tcal), tcalMean = mean(tcal), power = mean(vaccGood)),
          list(trial, gs, ord, delayUnit)]
## finInfo has 4 categories for each simulation in finTrials (all=all cases by trial end date,
## analyzed = all analyzed cases by trial end date, allFinalEV/no_EV = all cases by stop date w/ or
## w/o end vaccine rollout)

####################################################################################################
## extract counterfactuals
load('data/vaccProp1.Rdata')
vaccProp <- vaccProp1
vaccProp[, simNum:=1:length(vaccEff)]

thing <- 'Equip-cfs-3pit'
## fincfs <- extractCFs(thing, verb=0)
load(file=file.path('BigResults',paste0(thing, '.Rdata')))
fincfs <- merge(fincfs, vaccProp, by = 'simNum', all.y=F) ## copy vaccProp in there (should do this in analysisFuns.R later
class(fincfs$cc) <- 'numeric'

## Examine
fincfs[, length(vaccEff), list(simNum,cf)]
finTrials[, length(vaccEff), simNum] ## 4 types of simulations
finTrials[, list(numSim = length(vaccEff)), list(simNum, gs, ord, trial, propInTrial)] ## these are the groupings
finTrials[, list(numSim = length(vaccEff)), list(simNum, gs, ord, trial, propInTrial)][,max(numSim)]
finInfo[cat=='allFinalEV', length(caseTot), list(simNum, gs, ord, trial, propInTrial)] ## 4 types of simulations

setkey(finInfo, gs, ord, trial, propInTrial, simNum, cat)
setkey(finTrials, gs, ord, trial, propInTrial, simNum)

## Merge them so we have simulation level-info (speed, result, etc) in each finInfo category
finit <- finInfo[finTrials[, list(gs, ord, trial, propInTrial, simNum, tcal, vaccCases, contCases, vaccGood, vaccBad, vaccEff, PHU)]]
finit[, lab:=factor(paste0(trial, c('','gs-')[as.numeric(gs==T)+1],'-', ord))] ## make useful labels

## Merge finit with fincfs so we can compare counterfactuals and factuals
fall <- rbind(fincfs[, list(lab = cf, cat = 'allFinalEV', propInTrial, caseTot, vaccEff, simNum)],
              finit[, list(lab, cat, propInTrial, caseTot, vaccEff, simNum)]) ##cat %in% c('all','allFinalEV','allFinal_noEV')
fall[, posv:= vaccEff>0] ## positive vaccine efficacy simulations


####################################################################################################
## Figures
####################################################################################################
pdf('Figures/case dens.pdf')
ggplot(finit[cat=='allFinalEV'], aes(caseTot, colour = lab)) +
  geom_density() + facet_wrap(~propInTrial, ncol=1)
dev.off()

pdf('Figures/caseXtime.pdf')
p <- ggplot(finit[cat=='allFinalEV' & gs==T], aes(as.numeric(tcal), caseTot, colour = vaccGood))
p + geom_point(alpha = 1/5) + facet_wrap(~propInTrial+ord, ncol=2)
dev.off()

pdf('Figures/caseXtime.pdf')
p <- ggplot(finit[propInTrial==.2 & cat=='allFinalEV' & gs==T], aes(as.numeric(tcal), caseTot, colour = vaccGood))
p + geom_point(alpha = 1/5) + facet_grid(posv~ord)
dev.off()

pdf('Figures/caseXtime end vs not.pdf')
p <- ggplot(finit[propInTrial==.2 & cat %in% c('all','allFinalEV') & gs==T], aes(as.numeric(tcal), caseTot, colour = vaccGood))
p + geom_point(alpha = 1/5) + facet_grid(posv~ord + cat)
dev.off()

pdf('Figures/tp box.pdf')
p <- ggplot(finit[cat=='allFinalEV'], aes(lab, caseTot))
p + geom_boxplot() + facet_wrap(~propInTrial, ncol=1)
dev.off()

pdf('Figures/cfs.pdf')
ggplot(fall[!lab %in% c('RCT-none','RCT-TU')], aes(caseTot, colour = lab)) +
  geom_density(lwd=2) + facet_wrap(~propInTrial, ncol=1) + scale_x_continuous(limits = c(0, 400))
dev.off()

pdf('Figures/cfs pos vacceff.pdf')
ggplot(fall[ !lab %in% c('RCT-none','RCT-TU') & cat =='allFinalEV'], aes(caseTot, colour = lab)) +
  geom_density(lwd=2) + facet_grid(posv ~  propInTrial) + scale_x_continuous(limits = c(0, 400))
dev.off()

pdf('Figures/cfs pos vacceff end vs final.pdf')
ggplot(fall[ !grepl('pop',lab) & posv==T &  !lab %in% c('RCT-none','RCT-TU') ], aes(caseTot, colour = lab)) +
  geom_density(lwd=1) + facet_grid(cat ~  propInTrial) + scale_x_continuous(limits = c(0, 400))
dev.off()

pdf('Figures/cfs.pdf')
ggplot(fall[!lab %in% c('RCT-none','RCT-TU')], aes(caseTot, colour = lab)) +
  geom_density(lwd=2) + facet_wrap(~propInTrial, ncol=1) + scale_x_continuous(limits = c(0, 400))
dev.off()

finit[cat=='allFinalEV'][1:2]

####################################################################################################
## Difference between matched counterfactuals and factuals

## Merge (big merger but worth it in speed later)
fin <- merge(finit, fincfs, by ='simNum', all.y=F, allow.cartesian = T)
class(fin$simNum) <- class(fin$cc)<- 'integer'
setkey(fin, gs, ord.x, simNum, cf, cat)

## Get difference between caseTot
ctot <- fin[cat %in% c('allFinalEV', 'allFinal_noEV'), 
            list(cc, vaccGood, vaccBad, tcal, caseTotvsCF = caseTot.x - caseTot.y, caseTotF = caseTot.x, caseTotCF = caseTot.y),             list(gs, ord.x, simNum, cf, cat)]
setkey(ctot, simNum, cc, cf)
satr <- CJ(simNum = 1:max(ctot$simNum), cc = 1:max(ctot$cc), cf = c('NTpop','VRpop'))
ctot <- ctot[satr, allow.cartesian=T]## Merge to make sure that empty columns exist
ctot

## SEGAULT HERE
## Get difference between caseTot
## [1] below is because rows are replicated by cf, except for cf & caseTot.y & SAE
ctot <- fin[cat %in% c('allFinalEV', 'allFinal_noEV'), 
            list(cc[1], vaccGood[1], vaccBad[1], tcal[1],
                 caseTotF = caseTot.x[1],
                 caseTotNT = caseTot.y[cf=='NTpop'],
                 caseTotVR = caseTot.y[cf=='VRpop']), 
            list(gs, ord.x, simNum, cat)]
ctot

## ctot[, list(ctF = mean(caseTotF),
##             ctCF = mean(caseTotCF),
##             diff = mean(caseTotvsCF)),
##      list(gs, ord.x, simNum, cf, cat)]

## cs <- ctot[order(cat), list(ctF = mean(caseTotF), ctCF = mean(caseTotCF), diff = mean(caseTotvsCF), 
##                             power = mean(vaccGood), tcal = mean(tcal)), 
##            list(gs, ord.x, cf, cat)]
## save(cs, file = file.path('Results','cs.Rdata'))

## load(file = file.path('Results','cs.Rdata'))
## class(cs$gs) <- 'logical'
## setkey(cs, gs, ord.x)

## cs[, lab:=paste0(c('gs','')[as.numeric(gs) + 1], 'RCT-',c('none','TU')[as.numeric(ord.x=='TU') + 1])]
## cs[, lab:=factor(lab)]

## pdf('Figures/test.pdf')
## plot(0,0, type = 'n', bty = 'n', ylim = c(0,.5), xlim = cs[.(T, 'TU'), range(ctF,ctCF)], las = 1, xlab='cases', ylab = 'power')
## cs[cf=='NTpop' & cat=='allFinalEV', points(ctF, power, pch = 15, cex = 2, col = as.numeric(lab))]
## legend('bottom', leg = cs[, unique(lab)], col = cs[, as.numeric(unique(lab))], bty = 'n', pch = 15)
## abline(v = cs[cf=='NTpop', ctCF], lty=2)
## abline(v = cs[cf=='VRpop', ctCF], lty=2)
## dev.off()

## pdf('Figures/test.pdf')
## plot(0,0, type = 'n', bty = 'n', ylim = c(0,.5), xlim = cs[.(T, 'TU'), range(ctF,ctCF)], las = 1, xlab='cases', ylab = 'power')
## cs[cf=='NTpop' & cat=='allFinalEV', points(ctF, power, pch = 15, cex = 2, col = as.numeric(lab))]
## legend('bottom', leg = cs[, unique(lab)], col = cs[, as.numeric(unique(lab))], bty = 'n', pch = 15)
## abline(v = cs[cf=='NTpop', ctCF], lty=2)
## abline(v = cs[cf=='VRpop', ctCF], lty=2)
## dev.off()

## break down what's happening before trial & after trial (due to roll out)
## parse out false pos/negative from correct identifications
## is there a simple linear relationship between equipoise & power? (medium risk vs high risk group, do you lose more info than equipoise gained)
## post-vacc rollout period may modulate linearity of this relationship
