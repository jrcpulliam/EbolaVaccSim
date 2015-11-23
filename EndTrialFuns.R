
## Vaccinate control groups once vaccine efficacy identified. Wrapper around below functions
endT <- function(parms) {
    if(parms$verbose>30) browser()
    if(with(parms, delayUnit==0 & tail(intStats$vaccGood,1))) { ## if the vaccine is good on the stop date
        ## instantly vaccinate everyone 1 week after trial ends
        parms <- within(parms, {
            EVpopH <- copy(popH)
            EVpopH[vaccDay > endTrialDay, vaccDay := tail(intStats$tcal,1) + instVaccDelay]
        })
    }else{ ## run one of the below functions on them
        endFXN <- get(paste0('end', parms$trial))
        parms <- endFXN(parms) ## if vaccGood (if statements internal)
        parms <- badVaccEndT(parms) ## if vaccBad (if statements internal)
    }
    parms <- within(parms, {
        ## Reset immune indices
        EVpopH[, immuneDay := vaccDay + immunoDelay] 
        EVpopH[, vacc := day >= vaccDay]
        EVpopH[, immune := day >= immuneDay]
        ## copy pop, but change vaccDays & immuneDays for simInfection below
        EVpop <- copy(pop) 
        EVpop[, vaccDay := EVpopH[day==0, vaccDay]]
        EVpop[, immuneDay := EVpopH[day==0, immuneDay]]
    })
    ## do infection process again post-end of trial if not SWCT (which proceeds as normal)
    if(with(parms, tail(intStats[, vaccGood | vaccBad],1))) { 
        parms <- simInfection(parms, whichDo = 'EVpop', startInfectingDay = parms$endTrialDay)
    }
    ## calculate # infected for various permutations
    parms <- makeSurvDat(parms, whichDo='EVpop') ## don't make stActiveEV because not analyzing (just vacc rollout data)
    return(parms)
}

cfSims <- function(parms) {
    if(parms$doCFs) {
        ## *V*accine *r*ollout simulation; *N*o *t*rial counter-factual simulation
        parms <- within(parms, {
            NTpopH <- copy(popH)
            NTpopH <- within(NTpopH, {
                infectDay <- vaccDay <- immuneDay <- Inf
                vacc <- immune <- F
            })
            VRpop <- NTpop <- copy(pop)
            NTpop <- within(NTpop, {
                infectDay <- vaccDay <- immuneDay <- Inf
            })
        })
        ## ##############################
        ## risk-prioritized roll-out, i.e. no randomization
        ## force order to be that of TU-risk-prioritized-RCT
        parmsVR <- copy(parms)
        if(parmsVR$ord!='TU') {
            parmsVR$ord <- 'TU'
            parmsVR <- reordRCT(parmsVR)
        }
        parms$VRpopH <- parmsVR$popH ## to avoid rewriting reordFxns with whichDo
        rm(parmsVR)
        parms <- setSWCTvaccDays(parms, whichDo='VRpop') ## set vaccination days
        parms <- setImmuneDays(parms, whichDo='VRpop')   ## set immune days
        for(cf in c('NTpop','VRpop')) {
            parms <- simInfection(parms, whichDo = cf, startInfectingDay = 0)
            parms <- makeSurvDat(parms, whichDo=cf)
        }
    }
    return(parms)
}
## parms$VRpopH[idByClus %in% c(1,151)]

badVaccEndT <- function(parms) within(parms, {
    if(tail(intStats$vaccBad,1)) { ## Vaccine doesn't work
        EVpopH <- copy(popH)
        EVpopH[vaccDay > endTrialDay, vaccDay := Inf] ## no more vaccinating after the trial is over
    }
})

## SWCT is fastest possible vaccination already, so no change
endSWCT <- function(parms) within(parms, {
    EVpopH <- copy(popH)
})

endCRCT <- function(parms) within(parms, {
    if(verbose>33) browser()
    EVpopH <- copy(popH)
    if(tail(intStats$vaccGood,1)) {
        ## remaining clusters get vaccinated starting delayUnit after endTrialDay or when the last cluster 
        vaccClusters <- EVpopH[vaccDay <= endTrialDay  & vaccDay!=Inf , unique(cluster)]            
        notYetVaccClusters <- EVpopH[vaccDay > endTrialDay & vaccDay!=Inf , unique(cluster)]
        contClusters <- EVpopH[vaccDay==Inf , unique(cluster)]
        unVaccClusters <- c(notYetVaccClusters, contClusters)
        vaccDaysLeft <- daySeqLong[daySeqLong >= endTrialDay + delayUnit]
        if(ord == 'none') ## same random order as initially set
            EVpopH[vaccDay==Inf, vaccDay := vaccDaysLeft[1] + delayUnit*(cluster - numClus/2 - 1)] 
        if(ord == 'BL') { ## vaccinate remaining clusters in order of highest to lowest based on trial end time
            EVclusIncRank <- EVpopH[cluster %in% unVaccClusters & idByClus==1 & day == endTrialDay,
                                    rev(order(clusHaz))]
        }
        if(ord == 'TU') { ## vaccinate clusters in order of time-updated highest to lowest
            EVclusIncRank <- NULL
            for(ii in 1:length(unVaccClusters)) {
                dd <- vaccDaysLeft[ii]
                tmpRank <- popHearly[cluster %in% unVaccClusters & idByClus==1 & day == dd - reordLag,
                                  rev(order(clusHaz))]
                tmpRank <- tmpRank[!tmpRank %in% EVclusIncRank]
                EVclusIncRank <- c(EVclusIncRank, tmpRank[1])
                rm(tmpRank,dd)
            }
        }
        if(ord %in% c('BL','TU')) {
            for(ii in 1:length(EVclusIncRank))
                EVpopH[cluster==unVaccClusters[EVclusIncRank[ii]],
                       vaccDay := vaccDaysLeft[ii]]
            rm(ii, EVclusIncRank) }
    }
})



endFRCT <- endRCT <- function(parms) within(parms, {
    EVpopH <- copy(popH)
    if(tail(intStats$vaccGood,1)) {
        if(verbose>34) browser()
        ## remaining half get vaccinated in each cluster starting delayUnit day decision made to vaccinate controls
        vaccClusters <- EVpopH[vaccDay <= endTrialDay  & vaccDay!=Inf , unique(cluster)]            
        notYetVaccClusters <- EVpopH[vaccDay > endTrialDay & vaccDay!=Inf , unique(cluster)]

        ## ##################################################
        ## OPTION 1: Continue with cluster sequence, but vaccinate full clusters. Then go back and vacc partially vacc clusters.
        if(RCTendOption==1) {

            ## ##############################
            ## vaccinate everyone (not just half) in those clusters that haven't yet had anyone vaccinated moving forward
            lastPreassignedVaccDay <- popH[vaccDay!=Inf, max(vaccDay)]
            startVaccContDay <- max(lastPreassignedVaccDay, endTrialDay)
            vaccDaysLeft <- seq(min(daySeqLong[daySeqLong > startVaccContDay]), 700, by = delayUnit) ## remaining vaccination days (delayUnit apart)
            EVpopH[cluster %in% notYetVaccClusters & vaccDay==Inf,
                   vaccDay := delayUnit*(cluster-1)] ## same as those pre-assigned to be vaccinated 
            ## EVpopH[cluster %in% notYetVaccClusters & idByClus%in%c(1,151), list(cluster, vaccDay)]

            ## ##############################
            ## Vaccinate the rest of partially vaccinated clusters afterwards. First get order in which to do it, if necessary.
            if(ord == 'BL') ## reorder vaccination sequence by vaccinating highest hazard clusters as determined when trial ended
                EVclusIncRank <- EVpopH[cluster %in% vaccClusters & vaccDay == Inf & idByClus==1 & day == firstVaccDayAfterTrialEnd-reordLag,
                                        rev(order(clusHaz))]
            if(ord == 'TU') {  ## reorder vaccination sequence by vaccinating highest hazard partially vacc clusters as vaccination happens
                EVclusIncRank <- NULL
                for(ii in 1:length(vaccClusters)) {
                    dd <- vaccDaysLeft[ii]
                    tmpRank <- popHearly[cluster %in% vaccClusters & idByClus==1 & day == dd - reordLag,
                                         rev(order(clusHaz))]
                    tmpRank <- tmpRank[!tmpRank %in% EVclusIncRank]
                    EVclusIncRank <- c(EVclusIncRank, tmpRank[1])
                }
            }

            ## Then update the vaccination days based on the order
            if(ord == 'none')
                EVpopH[cluster %in% vaccClusters & vaccDay==Inf, 
                       vaccDay := vaccDaysLeft[1] + delayUnit*(cluster-min(vaccClusters))]
            if(ord %in% c('BL','TU')) { 
                for(ii in 1:length(EVclusIncRank))
                    EVpopH[cluster==EVclusIncRank[ii] & vaccDay==Inf, vaccDay := vaccDaysLeft[ii]]
                rm(ii) }
        }

        ## ##################################################
        ## Option 2. (Equipoise-maximizing) Vacc partially vacc clusters first, then vaccinate entire unvacc clusters.
        ## vaccinate the rest of partially vaccinated clusters afterwards
        if(RCTendOption==2) {
            ## Vaccinate the rest of partially vaccinated clusters first. 
            ## First get order in which to do it, if ord!='none'
            vaccDaysLeft <- seq(min(daySeqLong[daySeqLong > endTrialDay]), 700, by = delayUnit) ## remaining vaccination days (delayUnit apart)
            ## 700 ensures that everyone eventually gets vaccinated, even if after trial stops keeping track
            if(ord == 'BL') ## reorder vaccination sequence by vaccinating highest hazard clusters as determined when trial ended
                EVclusIncRank <- EVpopH[cluster %in% vaccClusters & idByClus==1 & day == firstVaccDayAfterTrialEnd-reordLag,
                                        rev(order(clusHaz))]
            if(ord == 'TU') {  ## reorder vaccination sequence by vaccinating highest hazard partially vacc clusters as vaccination happens
                EVclusIncRank <- NULL
                for(ii in 1:length(vaccClusters)) {
                    dd <- vaccDaysLeft[ii]
                    tmpRank <- popHearly[cluster %in% vaccClusters & idByClus==1 & day == dd - reordLag,
                                         rev(order(clusHaz))]
                    tmpRank <- tmpRank[!tmpRank %in% EVclusIncRank]
                    EVclusIncRank <- c(EVclusIncRank, tmpRank[1])
                }
            }
            ## Look at how this actually works
            ## popHearly[cluster %in% vaccClusters & idByClus==1 & day == dd - reordLag,
            ##                                       list(cluster, clusHaz, order(clusHaz), rev(order(clusHaz)))]

            ## ##############################
            ## Then update the vaccination days for unvaccinated individuals in partially vaccinated clusters based on these orders
            if(ord == 'none')
                EVpopH[cluster %in% vaccClusters & vaccDay==Inf, 
                       vaccDay := vaccDaysLeft[1] + delayUnit*(cluster-min(vaccClusters))]
            if(ord %in% c('BL', 'TU')) { 
                for(ii in 1:length(EVclusIncRank))
                    EVpopH[cluster==vaccClusters[EVclusIncRank[ii]] & vaccDay==Inf, vaccDay := vaccDaysLeft[ii]]
                rm(ii) }

            ## ##############################
            ## Then go back and vaccinate unvaccinated clusters
            if(length(notYetVaccClusters)>0) {
                ## Vacc days left after vaccinated control individuals in clusters that already had people vaccinted
                vaccDaysLeft <- seq(min(daySeqLong[daySeqLong > EVpopH[cluster %in% vaccClusters, max(vaccDay)]]), 700, by = delayUnit)
                ## Again get order first, for totally unvaccinated clusters
                if(ord == 'BL') ## reorder vaccination sequence by vaccinating highest hazard clusters as determined when *trial ended*
                    EVclusIncRank <- EVpopH[cluster %in% notYetVaccClusters & vaccDay == Inf & idByClus==1 & 
                                                day == firstVaccDayAfterTrialEnd-reordLag, rev(order(clusHaz))]
                if(ord == 'TU') { ## reorder vaccination sequence by vaccinating highest hazard remaining clusters as vaccination happens
                    EVclusIncRank <- NULL
                    for(ii in 1:length(notYetVaccClusters)) {
                        dd <- vaccDaysLeft[ii]
                        tmpRank <- popHearly[cluster %in% notYetVaccClusters & idByClus==1 & day == dd - reordLag,
                                             rev(order(clusHaz))]
                        tmpRank <- tmpRank[!tmpRank %in% EVclusIncRank]
                        EVclusIncRank <- c(EVclusIncRank, tmpRank[1])
                    }
                }
                ## Then set vaccination day for ALL within cluster (not those that
                ## haven't been already assigned, i.e. don't require vaccDay==Inf)
                if(ord == 'none')
                    EVpopH[cluster %in% notYetVaccClusters,
                           vaccDay := vaccDaysLeft[1] + delayUnit*(cluster-min(notYetVaccClusters))] 
                if(ord %in% c('BL','TU')) { 
                    for(ii in 1:length(EVclusIncRank))
                        EVpopH[cluster==notYetVaccClusters[EVclusIncRank[ii]],
                               vaccDay := vaccDaysLeft[ii]]
                    rm(ii) }
            }

            
        }

        ## ##################################################
        ## Option 3. (Incidence minimizing) Procede to vaccinate all unvaccinated individuals in order
        ## of their cluster hazard * (# participants lef unvaccinated)
        if(RCTendOption==3) {
            vaccDaysLeft <- seq(min(daySeqLong[daySeqLong > endTrialDay]), 700, by = delayUnit) ## remaining vaccination days 

            if(ord == 'BL') {## reorder vaccination sequence by vaccinating highest hazard clusters as determined when trial ended
                EVpopHclusTmp <- EVpopH[day == firstVaccDayAfterTrialEnd-reordLag, list(clusHaz = clusHaz[1]), cluster]
                EVpopHclusTmp[, stillUnInfUnVacc := 
                                  pop[, sum(infectDay > firstVaccDayAfterTrialEnd-reordLag & vaccDay > firstVaccDayAfterTrialEnd-reordLag),
                                      cluster]$V1]
                EVpopHclusTmp[, incRank := clusHaz * stillUnInfUnVacc]
                EVclusIncRank <- rev(order(EVpopHclusTmp$incRank))
            }

            if(ord == 'TU') {  ## reorder vaccination sequence by hazard*#unvacc at risk order
                EVclusIncRank <- NULL
                for(ii in 1:numClus) {
                    dd <- vaccDaysLeft[ii] - reordLag
                    EVpopHclusTmp <- EVpopH[day == dd, list(clusHaz = clusHaz[1]), cluster]
                    EVpopHclusTmp[, stillUnInfUnVacc := 
                                      pop[, sum(infectDay > dd & vaccDay > dd), cluster]$V1]
                    EVpopHclusTmp[, incRank := clusHaz * stillUnInfUnVacc]
                    tmpRank <- rev(order(EVpopHclusTmp$incRank))
                    tmpRank <- tmpRank[!tmpRank %in% EVclusIncRank]
                    EVclusIncRank <- c(EVclusIncRank, tmpRank[1])
                }
            }

            ## ##############################
            ## Proceed to finish vaccinating everyone based on their order
            if(ord == 'none')
                EVpopH[vaccDay==Inf, 
                       vaccDay := vaccDaysLeft[1] + delayUnit*(cluster-1)]
            if(ord %in% c('BL', 'TU')) { 
                for(ii in 1:length(EVclusIncRank))
                    EVpopH[cluster==c(1:numClus)[EVclusIncRank[ii]] & vaccDay==Inf, vaccDay := vaccDaysLeft[ii]]
                rm(ii) }

        }
    } ## End vaccGood if statement
    ## check ordering works
    ## EVpopH[cluster %in% notYetVaccClusters & idByClus %in% c(1,clusSize/2+1) & day == 0, list(cluster, idByClus, clusHaz, vaccDay)]
    ## EVpopH[cluster %in% vaccClusters & idByClus %in% c(1,clusSize/2+1) & day == 0, list(cluster, idByClus, clusHaz, vaccDay)]
})


## Create summary of # of cases by whether they were currently considered immune or not, & whether
## they were assigned to vaccination or control groups at baseline.
makeCaseSummary <- function(parms) within(parms, {
    if(parms$verbose>32) browser()
    stNms <- c('st','EVst','stActive','EVstActive')
    ## Create final survival tables for all of these censoring at trial's completion (not end day, but final time planned)
    ## ##################################################
    ## by person time spent in immune category (not omnitient immune but thought post vacc eclispe period)
    ptByPTImmune <- data.table(immuneGrp = 0:1)
    casesByPTImmune <- data.table(immuneGrp = 0:1)
    for(stI in stNms) {
        tmp <- censSurvDat(parms, whichDo = stI)
        casesTmp <- summarise(group_by(tmp, immuneGrp), sum(infected), sum(perstime))
        ptTmp <- arrange(casesTmp, immuneGrp)[,3, with=F]
        casesTmp <- arrange(casesTmp, immuneGrp)[,2, with=F]
        casesByPTImmune[,(stI) := casesTmp]
        nmTmp <- paste0(stI,'_pt')
        ptByPTImmune[,(nmTmp) := ptTmp]
    }
    casesXPT_Immune <- data.table(type = stNms) 
    casesXPT_Immune[, susCases := t(casesByPTImmune)[-1,1]]
    casesXPT_Immune[, immCases := t(casesByPTImmune)[-1,2]]
    casesXPT_Immune[, susPT := t(ptByPTImmune)[-1,1]]
    casesXPT_Immune[, immPT := t(ptByPTImmune)[-1,2]]
    ## ##################################################
    ## by randomization group
    ptByVaccRand <- data.table(vaccRandGrp = 0:1)
    casesByVaccRand <- data.table(vaccRandGrp = 0:1)
    vaccRandIndiv <- pop[vaccDay!=Inf, unique(indiv)]
    for(stI in stNms) {
        tmp <- censSurvDat(parms, whichDo = stI)
        tmp$vaccRandGrp <- tmp[, indiv %in% vaccRandIndiv]
        casesTmp <- summarise(group_by(tmp, vaccRandGrp), sum(infected), sum(perstime))
        ptTmp <- arrange(casesTmp, vaccRandGrp)[,3, with=F]
        casesTmp <- arrange(casesTmp, vaccRandGrp)[,2, with=F]
        casesByVaccRand[,(stI) := casesTmp]
        nmTmp <- paste0(stI,'_pt')
        ptByVaccRand[,(nmTmp) := ptTmp]
    }
    casesXVaccRandGrp <- data.table(type = stNms) 
    if(trial=='SWCT') { ## no one randomized to control group in SWCT
        casesXVaccRandGrp[, contCases := 0]
        casesXVaccRandGrp[, contPT := 0]
    }
    if(!trial=='SWCT') {
        casesXVaccRandGrp[, contCases := t(casesByVaccRand)[-1,1]]
        casesXVaccRandGrp[, contPT := t(ptByVaccRand)[-1,1]]
    }
    casesXVaccRandGrp[, vaccCases := t(casesByVaccRand)[-1,2]]
    casesXVaccRandGrp[, vaccPT := t(ptByVaccRand)[-1,2]]
    rm(nmTmp, tmp, casesTmp, ptTmp, stI, casesByPTImmune, ptByVaccRand,casesByVaccRand, ptByPTImmune, stNms, vaccRandIndiv)
})
