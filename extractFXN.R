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
    ## finTrials[,list(vaccEff,mean,PHU)] ## NEED TO AVERAGE BIAS ON PHU scale 

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
