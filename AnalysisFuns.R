testZeros <- function(parmsTmp) {
    tmpCSD <- censSurvDat(parmsTmp)
    casesXgroup <- tmpCSD[,list(cases = sum(infected)), immuneGrp]
    return(0 %in% casesXgroup[,cases])
}

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

prepDat <- function(parms, bump=T) {
    if(gs) { ## prepare group sequential analysis times
        parms <- within(parms, {
            maxInfo <- 35 ## number of events considered for sufficient power
            ## # events at which analyses occur
            intTab <- data.table(events = round(gsBounds$timing * maxInfo))
            ## Get vector of event (infection) timings
            infDays <- stActive$infectDay[stActive$infectDay!=Inf]
            infDays <- infDays[order(infDays)]
            ## Calculate interim analyses times
            intTab[, interimTimes:= ceiling(infDays[infoT])]
            intTab <- intTab[!is.na(interimTimes)]
            intTab$trigger <- 'events'
            ## If don't do all analyses, add one more at maximum trial duration
            if(nrow(intTab) < gsBounds$k) intTab <- rbind(intTab, data.table(events=NA, interimTimes=maxInfectDay, trigger='end time'))
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
            print(intTab)
            if(verbose==2.89) browser()
        })
    }
    
    for(tc in 1:length(timing)) {


    }
    if(!testZeros(parms)) { ## zero events in either arm?
        parmsE <- parms
        parmsE$bump <- F
    }else{
        parmsE <- infBump(parms)
        parmsE$bump <- T
    }
    tmpCSDE <- tmpCSD <- censSurvDat(parms)
    tmpASD <- censSurvDat(parms, whichDo='st') ## all survival data (not just actively analyzeable person-time), censored by trial end date
    if(testZeros(parms))  tmpCSDE <- censSurvDat(parmsE)

}

getEndResults <- function(parms, bump = T) {
    if(!testZeros(parms)) {
        parmsE <- parms
        parmsE$bump <- F
    }else{
        parmsE <- infBump(parms)
        parmsE$bump <- T
    }
    tmpCSDE <- tmpCSD <- censSurvDat(parms)
    tmpASD <- censSurvDat(parms, whichDo='st') ## all survival data (not just actively analyzeable person-time), censored by trial end date
    if(testZeros(parms))  tmpCSDE <- censSurvDat(parmsE)
    within(parmsE, {
        if(verbose==2.9) browser()
        vaccEE_ME_GS <- doCoxME_GS(parmsE, tmpCSDE, bump = bump)
        vaccEE_ME <- doCoxME(parmsE, tmpCSDE, bump = bump)
        ## vaccEE_GEEclusAR1 <- doGEEclusAR1(clusDat, csd=tmpCSDE, bump = bump)
        ## vaccEE_GLMMclus <- doGLMMclus(parmsE,, csd=tmpCSDE, bump = bump)
        ## vaccEE_GLMclus <- doGLMclus(parmsE, csd=tmpCSDE, bump = bump)
        vaccEE_GLMFclus <- doGLMFclus(parmsE, csd=tmpCSDE, bump = bump)
        ## vaccEE_MErelab <- doRelabel(parms, csd=tmpCSD, bump=F, nboot=nboot, verbFreqRelab=10)
        ## vaccEE_MEboot <- doBoot(parms, csd=tmpCSD, bump=F, nboot=nboot, verbFreqBoot=10)
        vEEs <- list(vaccEE_ME
                     ## , vaccEE_GLMMclus
                     ## , vaccEE_GLMclus
                   , vaccEE_GLMFclus
                     ## , vaccEE_GEEclusAR1
                     ## , vaccEE_MEboot
                     ## , vaccEE_MErelab
                     )
        finMods <- rbindlist(vEEs)
        finInfo <- compileStopInfo(tmp = tmpCSD, minDay=maxInfectDay,  verbose=verbose) ## active person-time only
        names(finInfo)[-1] <- paste0(names(finInfo)[-1],'_Active')
        finInfo <- data.frame(finInfo,  ## all person-time, not just active
                              compileStopInfo(tmp = tmpASD, minDay=maxInfectDay,  verbose=verbose)[-1])
        rm(vaccEE_ME
           ## , vaccEE_MEboot, vaccEE_MErelab
           ## , vaccEE_GEEclusAR1
           ## , vaccEE_GLMMclus 
         , vaccEE_GLMFclus
           ## , vaccEE_GLMclus
         , vEEs
           )
        return(list(finInfo=finInfo, finMods=finMods))
    })
}

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
        ## plotSTA(res$stActive) ## look at person-time for each data structure
        ## plotClusD(res$clusD)
        res <- prepDat(res)
        res <- getEndResults(res)
        finTmp <- data.frame(sim = ss, res$finMods)
        finMods <- rbind(finMods, finTmp)
        finITmp <- data.frame(sim = ss, res$finInfo)
        finInfo <- rbind(finInfo, finITmp)
        if(doSeqStops) {
            res <- seqStop(res)
            ## if(showSeqStops) {
            ##     resfull <- seqStop(res, fullSeq = T)
            ##     showSeqStop(resfull)
            ## }
            res <- endT(res)
            res <- makeCaseSummary(res)
            stopPt <- as.data.frame(tail(res$weeklyAns,1)) ## active cases by immmune grouping at time of case at end of trial
            stopPt <- with(res, {
                cbind(stopPt
                      , caseCXrandFinA = casesXVaccRandGrp[type=='EVstActive', contCases] ## active cases by vaccination randomization group at final
                      , caseVXrandFinA = casesXVaccRandGrp[type=='EVstActive', vaccCases]
                      , hazCXrandFinA = casesXVaccRandGrp[type=='EVstActive', contCases/contPT]/yearToDays
                      , hazVXrandFinA = casesXVaccRandGrp[type=='EVstActive', vaccCases/vaccPT]/yearToDays
                      , caseCXrandFin = casesXVaccRandGrp[type=='EVst', contCases]         ## total cases by vaccination randomization group at final
                      , caseVXrandFin = casesXVaccRandGrp[type=='EVst', vaccCases]
                      , hazCXrandFin = casesXVaccRandGrp[type=='EVst', contCases/contPT]/yearToDays
                      , hazVXrandFin = casesXVaccRandGrp[type=='EVst', vaccCases/vaccPT]/yearToDays
                      )
            })
            stopPoints <- rbind(stopPoints, stopPt)
        }
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


