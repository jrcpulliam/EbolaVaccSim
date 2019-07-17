if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
library(RColorBrewer); library(boot)
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

thing <- 'Equip1'
batchdirnm <- file.path('BigResults',thing)
fls <- list.files(batchdirnm, pattern='.Rdata', full.names = T)
length(fls)

dparms <- c('trial','gs','vaccEff','doSL','propInTrial','nbsize','ord','reordLag','delayUnit'##,'immunoDelay','trialStartDate'
            ##, 'weeklyDecay', 'cvWeeklyDecay', 'cvClus', 'cvClusTime',
            )
nbatch <- length(fls)
finInfoList <- finModList <- stopList <- parmsList <- list(NULL)
for(ii in 1:nbatch) {
    if(ii%%100 ==0) print(ii)
    ff <- fls[ii]
    if(exists('sim')) rm(sim)
    load(ff)
    if(exists('sim')) {
        sim$parms[['trialStartDate']] <- as.character(sim$parms[['trialStartDate']])
        parmsList[[ii]] <- data.frame(nbatch = ii, t(unlist(sim$parms[dparms])))
        tmpMod <- data.frame(nbatch = ii, sim$sim$finMods)
        finModList[[ii]] <- merge(tmpMod, sim$sim$finInfo, by = 'sim')
    }
}

parmsDT <- rbindlist(parmsList, use.names = T, fill = T)

finTrials <- merge(rbindlist(finModList), parmsDT, by = c('nbatch'))
finTrials[,vaccEff := levels(vaccEff)[vaccEff]]

finTrials[cat=='analyzed', list(mean(tcal), mean(vaccGood)), list(vaccEff, trial, gs, ord, delayUnit)]

finTrials[tcal >168 & cat=='analyzed']

finTrials[, sum(is.na(p)), mod]
finTrials[, sum(is.na(lci)), mod]
finTrials[, length(lci), list(propInTrial, mod)]
finTrials[mod=='coxME' & is.na(p), err:=1] ## sometimes cox returns NaNs, or partial NA's for certain values
finTrials$vaccEff <- as.numeric(finTrials$vaccEff)
