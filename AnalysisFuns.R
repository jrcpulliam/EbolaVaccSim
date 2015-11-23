compileStopInfo <- function(minDay, tmp, verbose=0) {
    if(verbose==4) browser()
    out <- data.table(stopDay=minDay
                    , caseCXimmGrpEnd = tmp[immuneGrp==0, sum(infected)]
                    , caseVXimmGrpEnd = tmp[immuneGrp==1, sum(infected)]
                    , hazCXimmGrpEnd = tmp[immuneGrp==0, sum(infected)/sum(perstime)]
                    , hazVXimmGrpEnd = tmp[immuneGrp==1, sum(infected)/sum(perstime)]
                    , ptRatioCVXimmGrpEnd = tmp[immuneGrp==0, sum(perstime)] / tmp[immuneGrp==1, sum(perstime)]
                      )
    out <- as.data.frame(out)
    return(out)
}



gsTimeCalc <- function(parms) {
    parms <- within(parms, {
        if(gs) { ## prepare group sequential analysis times
            if(verbose==2.89) browser()
            ## # events at which analyses occur
            intTab <- data.table(events = round(gsBounds$timing * maxInfo))
            ## Get vector of event (infection) timings
            infDays <- stActive$infectDay[stActive$infectDay!=Inf]
            infDays <- infDays[order(infDays)]
            ## Calculate interim analyses times
            intTab[, tcal:= ceiling(infDays[events])] ## calendar time (as opposed to information time)
            intTab <- intTab[!is.na(tcal)]
            intTab$trigger <- 'events'
            ## If don't do all analyses, add one more at maximum trial duration
            if(nrow(intTab) < gsBounds$k) intTab <- rbind(intTab, data.table(events=NA, tcal=maxDurationDay, trigger='end time'))
            intTab
            if(nrow(intTab) < gsBounds$k) { ## if don't have full # of analyses, must readjust design to spend all remaining alpha at maximum trial duration
                gsDesArgsAdj <- within(gsDesArgs, {
                    k <- nrow(intTab)
                    timing <- c(timing[k-1],1)
                })
                if(gsDesArgsAdj$k>1) { ## if any interims
                    gsBoundsAdj <- do.call(gsDesign, gsDesArgsAdj)
                    intTab <- cbind(intTab, upperZ = gsBoundsAdj$upper$bound, lowerZ = gsBoundsAdj$lower$bound)
                }else{ ## otherwise non-sequential
                    intTab <- data.table(events = NA, tcal = maxDurationDay, trigger = 'end time', upperZ = qnorm(.975), lowerZ = qnorm(.025) )
                }
            }else{ ## use original boundsQ
                gsBoundsAdj <- gsBounds
                intTab <- cbind(intTab, upperZ = gsBoundsAdj$upper$bound, lowerZ = gsBoundsAdj$lower$bound)
            }
            rm(maxInfo, infDays)
        }else{ ## non-sequential design
            intTab <- data.table(events = NA, tcal = maxDurationDay, trigger = 'end time', upperZ = qnorm(.975), lowerZ = qnorm(.025) )
        }
        intTab$contCases <- intTab$vaccCases <- intTab$obsZ <- as.numeric(NA)
        intTab$vaccGood <- intTab$vaccBad <- as.logical(NA)
        if(verbose>3) print(intTab)
    })
    return(parms)
}

testZeros <- function(tmpCSD) {
    casesXgroup <- tmpCSD[,list(cases = sum(infected)), immuneGrp]
    return(0 %in% casesXgroup[,cases])
}

getEndResults <- function(parms, bump = T) {
    if(verbose==2.93) browser()
    ## initialize
    trialStopped <- F
    analysisNum <- 0
    parms$intStats <- list()
    ## loop over sequential analyses (only do loop once for non-sequential analysis)
    while(!trialStopped) { 
        analysisNum <- analysisNum+1 ## iterate
        if(verbose>1) print(paste0('interim analysis ', analysisNum, ' of ', nrow(parms$intTab)))
        analysisDay <- parms$intTab[analysisNum, tcal]
        tmpCSDE <- tmpCSD <- censSurvDat(parms, censorDay = analysisDay)
        if(verbose>2) print(tmpCSD[, list(numInfected=sum(infected)), immuneGrp])
        ## Bump in case of 0-event arms
        if(!testZeros(tmpCSD)) { ## >0 events in each arm
            parmsE <- parms
            parmsE$bump <- F
        }else{ ## at least 1 arm has 0 events
            parmsE <- infBump(parms, censorDay=analysisDay)
            parmsE$bump <- T
            tmpCSDE <- censSurvDat(parmsE, censorDay = analysisDay)
            if(verbose>2) print(tmpCSDE[, list(numInfected=sum(infected)), immuneGrp])
        }
        parms <- within(parms, {
            ## Call analysis functions
            intStats[[analysisNum]] <- doStats(parmsE, tmpCSDE, analysisNum=analysisNum)
            ## Use negative z, since we think about crossing upper Z threshold as identifying
            ## positive vaccine, yet HR < 1 is equivaleynt.  use first StatsFxns item to determine stopping
            ## (usually CoxME), could vectorize this later but confusing to have different vaccination rollout
            ## strategies for one simulation due to different stopping times by different analyses
            intTab[analysisNum, obsZ:= - intStats[[analysisNum]][sf==StatsFxns[1], z]]
            intStats[[analysisNum]] <- cbind(intStats[[analysisNum]], analysis = analysisNum, numAnalyses = nrow(parms$intTab), intTab[analysisNum])
            ## Even for non-GS analyses, using beta/se = Z for significance testing, equivalent to p value from CoxPH
            intStats[[analysisNum]][, vaccGood :=  intTab[analysisNum, obsZ > upperZ] ]
            intStats[[analysisNum]][, vaccBad :=  intTab[analysisNum, obsZ < lowerZ] ]
            intStats[[analysisNum]][, contCases := tmpCSD[immuneGrp==0,sum(infected)]]
            intStats[[analysisNum]][, vaccCases := tmpCSD[immuneGrp==1,sum(infected)]]
        })
        earlyStop <- parms$intStats[[analysisNum]][, vaccGood | vaccBad]
        ## Determine whether trial stopped for boundary crossing or last analysis
        if(earlyStop | analysisNum==nrow(parms$intTab)) trialStopped <- T
    }
    tmpASD <- censSurvDat(parms, censorDay=analysisDay, whichDo='st') ## all survival data (not just actively analyzeable person-time), censored by trial end date
    parms <- within(parms, {
        ## Make intStats into one data table
        if(gsDesArgs$k>1) {
            intStats <- rbindlist(intStats)
        }else{
            intStats <- intStats[[1]]
        }
        if(tail(intStats[, vaccGood | vaccBad],1)) {
            endTrialDay <- tail(intStats$tcal, 1) ## when trial stopped (two-sided)
            firstVaccDayAfterTrialEnd <- min(daySeqLong[daySeqLong>endTrialDay])
        }else{
            endTrialDay <- maxDurationDay ## or maximum duration
        }
        ## Compile info on cases,  person-time & hazard at trial stopday
        finInfo <- compileStopInfo(tmp = tmpCSD, minDay=endTrialDay,  verbose=verbose) ## active person-time only
        names(finInfo)[-1] <- paste0(names(finInfo)[-1],'_Active')
        finInfo <- data.frame(finInfo,  ## all person-time, not just active
                              compileStopInfo(tmp = tmpASD, minDay=maxDurationDay,  verbose=verbose)[-1])

    })
    return(parms)
}

doStats <- function(parmsE, tmpCSDE, analysisNum=1) {
    with(parmsE, {
        if(verbose==2.94) browser()
        vEEs <- list()
        length(vEEs) <- length(StatsFxns)
        for(sf.ind in 1:length(StatsFxns)) {
            tempsf <- get(StatsFxns[sf.ind])
            argList <- list(parms=parmsE, csd=tmpCSDE, bump=parmsE$bump, nboot=parmsE$nboot)
            argList <- subsArgs(argList, tempsf)
            vEEs[[sf.ind]] <- do.call(tempsf, args = argList)
        }
        tmpStat <- rbindlist(vEEs)
        tmpStat$sf <- StatsFxns
        return(tmpStat)
    })
}
## StatsFxns <- c('doCoxMe','doGLMFclus','doGMMclus','doGLMclus','doRelabel','doBoot')

simNtrials <- function(seed = 1, parms=makeParms(), N = 2, returnAll = F,
                       doSeqStops = F, showSeqStops = F, flnm='test', verbFreq=10) {
    set.seed(seed)
    caseXVaccRandGrpList <- caseXPT_ImmuneList <- weeklyAnsList <- list()
    length(caseXVaccRandGrpList) <- length(caseXPT_ImmuneList) <- length(weeklyAnsList) <- N
    finInfo <- finMods <- stopPoints <- data.frame(NULL)
    for(ss in 1:N) {
        if(parms$verbose>0 & (ss %% verbFreq == 0)) print(paste('on',ss,'of',N))
        if(parms$verbose>.5 & (ss %% 1 == 0)) print(paste('on',ss,'of',N))
        if(parms$verbose==2) browser()
        res <- simTrial(parms)
        res <- makeSurvDat(res)
        res <- makeGEEDat(res)
        res <- activeFXN(res)
        res <- gsTimeCalc(res)
        ## plotSTA(res$stActive) ## look at person-time for each data structure
        ## plotClusD(res$clusD)
        res <- getEndResults(res)
        res$verbose <- 35
        res$endTrialDay <- 100
endT(res)
browser()


        finTmp <- data.frame(sim = ss, res$intStats[[length(res$intStats)]][1,])
        finMods <- rbind(finMods, finTmp)
        finITmp <- data.frame(sim = ss, res$finInfo)
        finInfo <- rbind(finInfo, finITmp)
        if(returnAll) {
            weeklyAnsList[[ss]] <- as.data.frame(res$weeklyAns)
            caseXVaccRandGrpList[[ss]] <- as.data.frame(res$casesXVaccRandGrp)
            caseXPT_ImmuneList[[ss]] <- as.data.frame(res$casesXPT_Immune)
        }
        rm(res)
        gc()
    }
    if(returnAll)
        return(list(
            stopPoints = stopPoints
            , weeklyAnsList = weeklyAnsList
            , caseXVaccRandGrpList = caseXVaccRandGrpList
            , caseXPT_ImmuneList = caseXPT_ImmuneList
            , finPoint=finPoint
            ))
    if(!returnAll)
        return(list(stopPoints=stopPoints, finMods=finMods, finInfo=finInfo))
}


