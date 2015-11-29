## Extract from simulations from multiple cores
extractSims <- function(thing, 
                        dparms = c('trial','gs','doSL','propInTrial','nbsize',
                            'ord','reordLag','delayUnit' ,'immunoDelay','trialStartDate'
                          , 'weeklyDecay', 'cvWeeklyDecay', 'cvClus', 'cvClusTime'
                                   )
                      , doSave=T
                      , verbose = 0
                        ) {
    if(verbose==1) browser()
    batchdirnm <- file.path('BigResults',thing)
    fls <- list.files(batchdirnm, pattern='.Rdata', full.names = T)
    nbatch <- length(fls)
    print(paste0('extracting from ', nbatch, ' files'))
    finInfoList <- finModList <- parmsList <- list(NULL)
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
    if(verbose==2) browser()
    parmsDT <- rbindlist(parmsList, use.names = T, fill = T)
    finTrials <- merge(rbindlist(finModList), parmsDT, by = c('nbatch'))
    ## print(finTrials[,list(tcalMean = mean(tcal), power = mean(vaccGood), falsePos = mean(vaccGood|vaccBad)), 
    ##                 list(vaccEff, trial, gs, ord, delayUnit)]
    
    ## Coverage
    finTrials[, cvr := lci < vaccEff & uci > vaccEff]
    ## Bias, must be done on RH/(RH+1) scale to deal with Inf & 0 RH's
    finTrials$RH <- finTrials[, 1-mean]
    finTrials$PHU <- finTrials[, RH/(RH+1)] ## proportion of hazard unavoidable even with vaccination
    finTrials[RH==Inf, PHU:=1] ## otheriwse gives NaN for Inf/(Inf+1)
    finTrials[, stopped:=vaccGood|vaccBad]

    ## finTrials[,list(vaccEff,mean,PHU)] ## NEED TO AVERAGE BIAS ON PHU scale 
    print("Distribution of finTrials$err")
    print(finTrials[, table(err)]) ## -1 couldn't fit coxme, so fit coxph instead, 1 couldn't fit anything at all

    ## Reorder columns
    front <- c('mod','vaccGood','vaccBad')
    setcolorder(finTrials, c(front, setdiff(names(finTrials), front)))
    back <- c('nbatch','sim')
    setcolorder(finTrials, c(setdiff(names(finTrials), back), back))
    if(doSave) {
        save(finTrials, file=file.path('BigResults', paste0(thing, '.Rdata')))
        save(finTrials, file=file.path('Results', paste0(thing, '.Rdata')))
    }
    return(finTrials)
}

makePowTable <- function(finTrials, verbose=0, doSave=T) {
    if(verbose==1) browser()
    pf <- finTrials[, list(
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


    if(verbose==2) browser()

    coerce.numeric <- function(dat = pf, variables = c('propInTrial', 'delayUnit','immunoDelay')) within(dat, {
        for(vv in variables) {
            if(dat[,class(get(vv))]=='character') dat[,vv:= as.numeric(get(vv)), with=F]
            if(dat[,class(get(vv))]=='factor') dat[,vv:= as.numeric(levels(get(vv))[get(vv)]), with = F]
        }
        rm(vv)
    })
    coerce.numeric(pf)

    pf[,trial:=factor(trial)]
    ## Get bias from means done on a PHU scale
    pf$RHxPHUNAR <- pf[, -PHUNAR/(PHUNAR-1)]
    pf$meanXPHUNAR <- pf[, 1 - RHxPHUNAR]
    pf$biasNAR <- pf[, meanXPHUNAR - as.numeric(vaccEff)]

    ## Formatting stuff
    front <- c('mod','vaccEff','stoppedNAR','vaccGoodNAR','cvrNAR','biasNAR',
               'nsim','meanErr','propInTrial','vaccBad','cvr','stopped','vaccGood')
    setcolorder(pf, c(front, setdiff(names(pf), front)))
    pf <- pf[!(trial=='FRCT' & delayUnit==0) & !(ord=='TU' & delayUnit==0)] ## redundant
    pf$trialStartDate <- as.Date(pf$trialStartDate)
    pf[mod=='coxME', mod:='CoxME']
    pf$mod <- factor(pf$mod, levels=unique(pf$mod))
    pf$order <- pf$ord
    pf[delayUnit==0, order:='simultaneous instant']
    pf$design <- pf$trial
    levels(pf$design)[levels(pf$design) == 'SWCT'] <- 'SWT'
    levels(pf$order)[2] <- 'time-updated'

    pf[, pit:=factor(paste0(propInTrial*100,'%'))]
    pf[, pit:=factor(pit, levels = c('2.5%','5%','7.5%','10%'), ordered = T)]
    baseMods <- c('Cox PH Frailty'
                , 'Poisson GLM\n no cluster effects'
                , 'Poisson GLM \nwith fixed effects by cluster')
    pf$model <- pf$mod
    levels(pf$model) <- paste0(rep(c('', 'bootstrap over\n', 'permutation test over\n'),each=3), rep(baseMods,3))
    if(doSave)
        save(pf, file=file.path('Results',paste0('pf_',thing,'.Rdata')))
    return(pf)
}
