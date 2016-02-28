
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
            intTab[, tcal:= ceiling(infDays[events])] ## calendar time (as opposed to information time); *CEILING MEANS MIGHT HAVE MORE INFO*
            intTab <- intTab[!is.na(tcal)]
            intTab$trigger <- 'events'
            intTab <- intTab[tcal <= maxDurationDay]
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
    if(parms$doTrial) { ## if doing a trial
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
                ## 
                ## for GS analyses, using beta/se = Z for significance testing, equivalent to p value from CoxPH
                if(gs) {
                    intTab[analysisNum, obsZ:= - intStats[[analysisNum]][sf==StatsFxns[1], z]]
                    intStats[[analysisNum]] <- cbind(intStats[[analysisNum]], analysis = analysisNum, numAnalyses = nrow(parms$intTab), 
                                                     intTab[analysisNum])
                    intStats[[analysisNum]][, vaccGood :=  intTab[analysisNum, obsZ > upperZ] ]
                    intStats[[analysisNum]][, vaccBad :=  intTab[analysisNum, obsZ < lowerZ] ]
                }else{
                    intStats[[analysisNum]] <- cbind(intStats[[analysisNum]], analysis = analysisNum, numAnalyses = nrow(parms$intTab), 
                                                     intTab[analysisNum])
                    intStats[[analysisNum]][, vaccGood :=  intTab[analysisNum, p < .05 & mean>0] ]
                    intStats[[analysisNum]][, vaccBad :=  intTab[analysisNum, p < .05 & mean <0] ]
                }
                intStats[[analysisNum]][, contCases := tmpCSD[immuneGrp==0,sum(infected)]]
                intStats[[analysisNum]][, vaccCases := tmpCSD[immuneGrp==1,sum(infected)]]
            })
            earlyStop <- parms$intStats[[analysisNum]][, vaccGood | vaccBad]
            ## Determine whether trial stopped for boundary crossing or last analysis
            if(earlyStop | analysisNum==nrow(parms$intTab)) trialStopped <- T
        }
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
        })
    }else{
        parms$endTrialDay <- NA ## no trial so no endday
    }
    return(parms)
}

## tmp <- do.call(doRelabel, args=list(parmsE, tmpCSDE, parmsE$bump, nboot = 5))

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

compileStopInfo <- function(atDay, tmp, verbose=0) {
    if(verbose==4) browser()
    out <- data.table(atDay=atDay
                    , caseC = tmp[immuneGrp==0, sum(infected)]
                    , caseV = tmp[immuneGrp==1, sum(infected)]
                    , hazC = tmp[immuneGrp==0, sum(infected)/sum(perstime)]
                    , hazV = tmp[immuneGrp==1, sum(infected)/sum(perstime)]
                    , ptRatioCV = tmp[immuneGrp==0, sum(perstime)] / tmp[immuneGrp==1, sum(perstime)]
                    , nsae = tmp[,sum(SAE)]
                      )
    ## To check that SAE's are working:
    ## tmp[indiv%in%tmp[SAE>0, unique(indiv)], list(indiv, startDay,endDay,vaccDay,immuneDayThink, immuneGrp, SAE)][order(indiv)]
    out <- as.data.frame(out)
    return(out)
}

finInfoFxn <- function(parms) {
    tempFXN <- function(atDay, whichDo, verbose=parms$verbose)
        compileStopInfo(tmp=censSurvDat(parms, censorDay=atDay, whichDo=whichDo), 
                        atDay=atDay, verbose=verbose)
    if(parms$doTrial) {
        compTab <- data.table(atDay = with(parms, c(endTrialDay, trackUntilDay))[c(1,1,2,2)]
                            , whichDo = c('stActive', 'st','stEV', 'st')
                            , lab = c('analyzed','all','allFinalEV','allFinal_noEV')
                              )
    }else{
        compTab <- data.table(atDay = parms$trackUntilDay, whichDo = 'st', lab = 'allFinalEV')
    }
    for(ii in 1:nrow(compTab)) {
        finInfoTmp <- do.call(tempFXN, args = as.list(compTab[ii, list(atDay, whichDo)]))
        if(ii==1) finInfo <- finInfoTmp else finInfo <- rbind(finInfo, finInfoTmp)
    }
    finInfo$cat <- compTab$lab
    finInfo <- as.data.table(finInfo)[order(atDay)]
    ## vaccEff estimate should roughly equate to 
    ## 1-finInfo[1,hazV/hazC]
    finInfo$caseTot <- finInfo[,caseC+caseV]
    setcolorder(finInfo, c('cat','atDay','caseTot','caseC','caseV','nsae','hazC','hazV','ptRatioCV'))
    parms$finInfo <- finInfo
    return(parms)
}

indivLevRes <- function(parms) within(parms, {
    SpopH <- popH[idByClus==1,list(Oc, day, Date, clusHaz)] ## cluster chars
    toDo <- c('pop','EVpop')
    if(!doTrial)   toDo <- 'pop'
    for(tmpnm in toDo) {
            tmp <- get(tmpnm)
            tmpH <- get(paste0(tmpnm,'H'))
            ## remember hazards are daily so multiply time unit
            tmpH[,indivHazVacc := indivHaz * c(rep(1, sum(!immune)), rep(1-vaccEff, sum(immune))), indiv]
            tmp <- merge(tmp, tmpH[, list(cumHaz = sum(indivHazVacc*hazIntUnit)), indiv], by = 'indiv')
            tmpSM <- tmp[,list(indiv, cluster, Oi, Oc, indivRR, vaccDay, infectDay, SAE, cumHaz)] ## indiv events
            assign(paste0('S',tmpnm), tmpSM)
        }
    if(doTrial) {## merge EV with no_EV
        Spop <- merge(Spop, SEVpop[,list(indiv, vaccDay, infectDay, SAE, cumHaz)], suffixes = c('','_EV'), by = c('indiv'))
        rm(SEVpop)
    }
    rm(tmp, tmpH, tmpSM)
})

simNtrials <- function(batch = 1, parms=makeParms(), N = 2, 
                       returnEventTimes=T, simNums = ((batch-1)*nsims + 1):((batch-1)*nsims + N),
                       verbFreq=10, vaccProp=NA) {
    finInfo <- finMods <- data.frame(NULL)
    SpopH <- Spop <- data.table()
    for(ss in 1:N) {
        simNum <- simNums[ss]
        set.seed(simNum)
        if(parms$verbose>0 & (ss %% verbFreq == 0)) print(paste('on',ss,'of',N))
        if(parms$verbose>.5 & (ss %% 1 == 0)) print(paste('on',ss,'of',N))
        if(parms$verbose==2) browser()
        if(!is.na(vaccProp)[1]) { ## set vaccine properties to value from pre-determined bayesian prior deviate
            parms$vaccEff <- vaccProp[simNum, vaccEff]
            parms$pSAE <- vaccProp[simNum, pSAE]
        }
        res <- simTrial(parms)
        res <- makeSurvDat(res)
        res <- makeGEEDat(res) ## probably unnecessary to do these lines for !doTrial if tweak the code a bit later
        res <- activeFXN(res)
        res <- gsTimeCalc(res)
        ## plotSTA(res$stActive) ## look at person-time for each data structure
        ## plotClusD(res$clusD)
        res <- getEndResults(res)
        res <- endT(res)
        res <- finInfoFxn(res)
        res <- indivLevRes(res)
        ## compile results from the final interim analysis (or all statistical analyses for a single fixed design)
        finTmp <- data.table(sim = ss, simNum = simNum, 
                             vaccEff=parms$vaccEff, pSAE=parms$pSAE)
        if(res$doTrial) finTmp <- data.table(finTmp, res$intStats[analysis==max(res$intStats$analysis)]) ## add analysis if doing a study
        finMods <- rbind(finMods, finTmp) ## model table
        finITmp <- data.table(sim = ss, simNum = simNum, res$finInfo)
        finInfo <- rbind(finInfo, finITmp) ## cases table
        if(returnEventTimes) {
            Spop <- rbind(Spop, data.table(sim=ss, simNum=simNum, res$Spop))
            SpopH <- rbind(SpopH, data.table(sim=ss, simNum=simNum, res$SpopH))
        }
        ## res <- equiCalc(res)
        rm(res)
        gc()
    }
    Spops <- list(Spop=Spop,SpopH=SpopH)
    return(list(finMods=finMods, finInfo=finInfo, Spops = Spops))
}


## ####################################################################################################
## ## To show that RNGs line up between counterfactuals & factuals for first simulation & with
## ## vaccEff=0 Otherwise, individuals' infection patterns differ necessarily. Even with a single
## ## piece-wise exponential RNG (rather than the current for loop across pieces), we couldn't maintain
## ## individuals having the same infection status up until their hazard changes between
## ## counterfactuals because waiting times for the same seeds would differ for different piece-wise
## ## exponential hazard functions even if the first few pieces were identical & the waiting time was
## ## likely in those first few pieces.
## parms$vaccEff <- 0
## system.time(simCF <- simN_CFs(batch=batch, parms=parms, N = nsims, returnInfTimes = T, verbFreq=10, vaccProp=NULL))
## system.time(sim <- simNtrials(batch=batch, parms=parms, N = nsims, verbFreq=10, vaccProp=NULL))
## sim$finInfo[cat=='allFinal_noEV']
## simCF$EventTimes[order(cc,simNum), list(caseTot = sum(infectDay<Inf), saeTot = sum(SAE)), list(simNum, ss, batch, cc, cf)]
## ####################################################################################################
