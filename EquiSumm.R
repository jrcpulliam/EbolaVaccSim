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

dparms <- c('trial','gs','vaccEff','doSL','propInTrial','nbsize','ord','reordLag','delayUnit' ,'immunoDelay','trialStartDate'
            , 'weeklyDecay', 'cvWeeklyDecay', 'cvClus', 'cvClusTime'
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
        parmsList[[ii]] <- data.table(nbatch = ii, t(unlist(sim$parms[dparms])))
        finModList[[ii]] <- data.table(nbatch = ii, sim$sim$finMods)
        finInfoList[[ii]] <- data.table(nbatch = ii, sim$sim$finInfo)
    }
}

parmsDT <- rbindlist(parmsList, use.names = T, fill = T)

finTrials <- merge(rbindlist(finModList), parmsDT, by = c('nbatch'))
finTrials[,list(tcalMean = mean(tcal), power = mean(vaccGood), falsePos = mean(vaccGood|vaccBad)), list(vaccEff, trial, gs, ord, delayUnit)]
    
## Coverage
finTrials[, cvr := lci < vaccEff & uci > vaccEff]

## Bias, must be done on RH/(RH+1) scale to deal with Inf & 0 RH's
finTrials$RH <- finTrials[, 1-mean]
finTrials$PHU <- finTrials[, RH/(RH+1)] ## proportion of hazard unavoidable even with vaccination
finTrials[RH==Inf, PHU:=1] ## otheriwse gives NaN for Inf/(Inf+1)
finTrials[,list(vaccEff,mean,PHU)] ## NEED TO AVERAGE BIAS ON PHU scale 

## Reorder columns
front <- c('mod','vaccGood','vaccBad')
setcolorder(finTrials, c(front, setdiff(names(finTrials), front)))
back <- c('nbatch','sim')
setcolorder(finTrials, c(setdiff(names(finTrials), back), back))

save(finTrials, file=file.path('BigResults', paste0(thing, '.Rdata')))
save(finTrials, file=file.path('Results', paste0(thing, '.Rdata')))

load(file=file.path('BigResults',paste0(thing, '.Rdata')))

finTrials[, table(err)] ## 18 times couldn't fit coxme, so fit coxph instead
finTrials[order(gs), list(tcalMean = mean(tcal), power = mean(vaccGood)), list(vaccEff, trial, gs, ord, delayUnit)]
finTrials[tcal<168, list(nsim=length(tcal), tcalMean = mean(tcal), power = mean(vaccGood)), list(vaccEff, trial, gs, ord, delayUnit)]
finTrials[, stopped:=vaccGood|vaccBad]

powFin <- finTrials[, list(
                     nsim = length(stopped)
                    , stopped = mean(stopped)
                    , vaccGood = mean(vaccGood)
                    , stopDay = mean(tcal)
                    , caseTot = mean(vaccCases+contCases)
                    , caseC = mean(contCases)
                    , caseV = mean(vaccCases)
                    , cvr = mean(cvr)
                    , cvrNAR = mean(cvr, na.rm=T)
                    , mean = mean(mean)
                    , meanNAR = mean(mean, na.rm=T)
                    , vaccBad = mean(vaccBad)
                    , stoppedNAR = mean(stopped,na.rm=T)
                    , vaccGoodNAR = mean(vaccGood,na.rm=T)
                    , vaccBadNAR = mean(vaccBad,na.rm=T)
                    , PHUNAR = mean(PHU, na.rm=T)
                    , meanErr = mean(err)
                    , meanBump = mean(bump)
                    ## , caseTot = mean(caseCXimmGrpEnd + caseVXimmGrpEnd )
                    ## , caseC = mean(caseCXimmGrpEnd)
                    ## , caseV = mean(caseVXimmGrpEnd)
                    ),
                    list(vaccEff, trial, gs, propInTrial, ord, delayUnit, mod, immunoDelay,trialStartDate
                       , weeklyDecay, cvWeeklyDecay, cvClus) #, cvClusTime)
                    ]


powFin[,propInTrial:= as.numeric(levels(propInTrial)[propInTrial])]
powFin[,delayUnit:= as.numeric(levels(delayUnit)[delayUnit])]
powFin[,trial:=factor(trial)]
## Get bias from means done on a PHU scale
powFin$RHxPHUNAR <- powFin[, -PHUNAR/(PHUNAR-1)]
powFin$meanXPHUNAR <- powFin[, 1 - RHxPHUNAR]
powFin$biasNAR <- powFin[, meanXPHUNAR - as.numeric(vaccEff)]
powFin[vaccEff==.7, list(gs, meanNAR,meanXPHUNAR,vaccEff, biasNAR, cvrNAR)]


powFin


## Formatting stuff
front <- c('mod','vaccEff','stoppedNAR','vaccGoodNAR','cvrNAR','biasNAR',
'nsim','meanErr','propInTrial','vaccBad','cvr','stopped','vaccGood')
setcolorder(powFin, c(front, setdiff(names(powFin), front)))
pf <- data.table(powFin)
pf <- pf[!(trial=='FRCT' & delayUnit==0) & !(ord=='TU' & delayUnit==0)] ## redundant
pf$trialStartDate <- as.Date(pf$trialStartDate)
pf[mod=='coxME', mod:='CoxME']
pf$mod <- factor(pf$mod, levels=unique(pf$mod))
pf$order <- pf$ord
pf[delayUnit==0, order:='simultaneous instant']
pf$design <- pf$trial
levels(pf$design)[levels(pf$design) == 'SWCT'] <- 'SWT'
levels(pf$order)[2] <- 'time-updated'
pf[, immunoDelay:=as.numeric(levels(immunoDelay)[immunoDelay])]
pf[, pit:=factor(paste0(propInTrial*100,'%'))]
pf[, pit:=factor(pit, levels = c('2.5%','5%','7.5%','10%'), ordered = T)]
baseMods <- c('Cox PH Frailty'
              , 'Poisson GLM\n no cluster effects'
              , 'Poisson GLM \nwith fixed effects by cluster')
pf$model <- pf$mod
levels(pf$model) <- paste0(rep(c('', 'bootstrap over\n', 'permutation test over\n'),each=3), rep(baseMods,3))

save(pf, file=file.path('Results',paste0('powFin_',thing,'.Rdata')))

