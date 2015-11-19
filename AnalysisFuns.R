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
            maxInfo <- 35 ## number of events considered for sufficient power
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
            if(nrow(intTab) < gsBounds$k) intTab <- rbind(intTab, data.table(events=NA, tcal=maxInfectDay, trigger='end time'))
            intTab
            if(nrow(intTab) < gsBounds$k) { ## if don't have full # of analyses, must readjust design to spend all remaining alpha at maximum trial duration
                gsDesArgsAdj <- within(gsDesArgs, {
                    k <- numInterims
                    timing <- c(timing[k-1],1)
                })
                gsBoundsAdj <- do.call(gsDesign, gsDesArgsAdj)
                intTab <- cbind(intTab, upperZ = gsBoundsAdj$upper$bound, lowerZ = gsBoundsAdj$lower$bound)
            }else{
                gsBoundsAdj <- gsBounds
                intTab <- cbind(intTab, upperZ = gsBoundsAdj$upper$bound, lowerZ = gsBoundsAdj$lower$bound)
            }
        rm(maxInfo, infDays)
        }else{ ## non-sequential design
            intTab <- data.table(events = NA, tcal = maxInfectDay, trigger = 'end time', upperZ = qnorm(.975), lowerZ = qnorm(.025) )
        }
        intTab$obsZ <- as.numeric(NA)
        print(intTab)
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
            vaccGood <-  intTab[analysisNum, obsZ > upperZ] 
            vaccBad <-  intTab[analysisNum, obsZ < lowerZ] 
        })
        ## Determine whether trial stopped for boundary crossing or last analysis
        if(parms$vaccGood|parms$vaccBad | analysisNum==nrow(parms$intTab)) trialStopped <- T
    }
    tmpASD <- censSurvDat(parms, censorDay=analysisDay, whichDo='st') ## all survival data (not just actively analyzeable person-time), censored by trial end date
    parms <- within(parms, {
        finInfo <- compileStopInfo(tmp = tmpCSD, minDay=maxInfectDay,  verbose=verbose) ## active person-time only
        names(finInfo)[-1] <- paste0(names(finInfo)[-1],'_Active')
        finInfo <- data.frame(finInfo,  ## all person-time, not just active
                              compileStopInfo(tmp = tmpASD, minDay=maxInfectDay,  verbose=verbose)[-1])
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
browser()
        finTmp <- data.frame(sim = ss, res$finMods)
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


