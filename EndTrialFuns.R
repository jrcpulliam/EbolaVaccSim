
## Vaccinate control groups once vaccine efficacy identified. Wrapper around below functions
endT <- function(parms, browse=F) {
    if(browse) browser()
    if(with(parms, delayUnit==0 & vaccEffEst['p'] < .05 & vaccEffEst['lci']>0)) { 
        ## instantly vaccinate everyone 1 week after trial ends
        parms <- within(parms, {
            EVpopH <- copy(popH)
            EVpopH[vaccDay > endTrialDay, vaccDay := endTrialDay + instVaccDelay]
        })
    }else{ ## run one of the below functions on them
        endFXN <- get(paste0('end', parms$trial))
        parms <- endFXN(parms)
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
    if(with(parms, trial != 'SWCT' & vaccEffEst['p'] < .05 & vaccEffEst['lci']>0)) { 
        parms <- simInfection(parms, whichDo = 'EVpop', startInfectingDay = parms$endTrialDay)
    }
    ## calculate # infected for various permutations
    parms <- makeSurvDat(parms, whichDo='EVpop')
    parms <- activeFXN(parms, whichDo = 'EVst')
    return(parms)
}

##         pop[infectDay!=Inf, list(indiv, infectDay)]
##         arrange(popH[infectDay!=Inf, list(indiv, infectDay)], indiv)
## with(parms, {
## browser()

##         EVpop[infectDay!=Inf, list(indiv, infectDay)]
##         arrange(EVpopH[infectDay!=Inf, list(indiv, infectDay)], indiv)

## nms <- colnames(pop)
## tst <- setcolorder(copy(popH[infectDay!=Inf])[,nms, with=F], nms)
## identical(pop[,1,with=F], tst[,1,with=F])
## identical(pop,tst)
## pop[which(pop[,infectDay]!=tst[,infectDay])]
## })

## identical(pop, EVpop)
## nms <- colnames(pop)
## tst <- setcolorder(copy(popH[day==0])[,nms, with=F], nms)
## identical(pop[,1,with=F], tst[,1,with=F])
## identical(pop,tst)
## pop[which(pop[,infectDay]!=tst[,infectDay])]

## for(ii in 1:ncol(pop)) print(setdiff(pop[,ii,with=F], tst[,ii,with=F]))

## identical(popH, EVpopH)

## tmp <- censSurvDat(parms, 259,'EVstActive')
## tmp[immuneGrp==1, sum(infected)]
## tmp[immuneGrp==0, sum(infected)]

## Create summary of # of cases by whether they were currently considered immune or not, & whether
## they were assigned to vaccination or control groups at baseline.
makeCaseSummary <- function(parms, browse=F) within(parms, {
    if(browse) browser()
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

## SWCT is fastest possible vaccination already, so no change
endSWCT <- function(parms) within(parms, {
    EVpopH <- copy(popH)
})

endCRCT <- function(parms) within(parms, {
    EVpopH <- copy(popH)
    if(vaccEffEst['p'] <.05 & vaccEffEst['lci']>0) {
        ## remaining clusters get vaccinated starting delayUnit after endTrialDay or when the last cluster 
        vaccClusters <- EVpopH[vaccDay <= endTrialDay  & vaccDay!=Inf , unique(cluster)]            
        notYetVaccClusters <- EVpopH[vaccDay > endTrialDay & vaccDay!=Inf , unique(cluster)]
        contClusters <- EVpopH[vaccDay==Inf , unique(cluster)]
        unVaccClusters <- c(notYetVaccClusters, contClusters)
        vaccDaysLeft <- daySeq[daySeq >= endTrialDay + delayUnit]
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
                tmpRank <- EVpopH[cluster %in% unVaccClusters & idByClus==1 & day == dd,
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



endFRCT <- endRCT <- function(parms, browse=F) within(parms, {
    EVpopH <- copy(popH)
    if(vaccEffEst['p'] <.05 & vaccEffEst['lci']>0) {
        if(browse) browser()
        ## remaining half get vaccinated in each cluster starting delayUnit day decision made to vaccinate controls
        vaccClusters <- EVpopH[vaccDay <= endTrialDay  & vaccDay!=Inf , unique(cluster)]            
        notYetVaccClusters <- EVpopH[vaccDay > endTrialDay & vaccDay!=Inf , unique(cluster)]

        ## ##################################################
        ## OPTION 1: Continue with cluster sequence, but vaccinate full clusters. Then go back and vacc partially vacc clusters.
        if(RCTendOption==1) {
            ## vaccinate everyone (not just half) in those clusters that haven't yet had anyone vaccinated moving forward
            lastPreassignedVaccDay <- popH[vaccDay!=Inf, max(vaccDay)]
            startVaccContDay <- max(lastPreassignedVaccDay, endTrialDay)
            vaccDaysLeft <- seq(min(daySeq[daySeq > startVaccContDay]), 700, by = delayUnit) ## remaining vaccination days (delayUnit apart)
            EVpopH[cluster %in% notYetVaccClusters & vaccDay==Inf,
                   vaccDay := delayUnit*(cluster-1)] ## same as those pre-assigned to be vaccinated
            ## Vaccinate the rest of partially vaccinated clusters afterwards. First get order in which to do it, if necessary.
            if(ord == 'BL') ## reorder vaccination sequence by vaccinating highest hazard clusters as determined when trial ended
                EVclusIncRank <- EVpopH[cluster %in% vaccClusters & vaccDay == Inf & idByClus==1 & day == endTrialDay,
                                        rev(order(clusHaz))]
            if(ord == 'TU') {  ## reorder vaccination sequence by vaccinating highest hazard partially vacc clusters as vaccination happens
                EVclusIncRank <- NULL
                for(ii in 1:length(vaccClusters)) {
                    dd <- vaccDaysLeft[ii]
                    tmpRank <- EVpopH[cluster %in% vaccClusters & idByClus==1 & day == dd,
                                      rev(order(clusHaz))]
                    tmpRank <- tmpRank[!tmpRank %in% EVclusIncRank]
                    EVclusIncRank <- c(EVclusIncRank, tmpRank[1])
                }
            }
            ## Then update the vaccination days based on the order
            if(ord == 'none')
                EVpopH[cluster %in% vaccClusters & vaccDay==Inf, 
                       vaccDay := vaccDaysLeft[1] + delayUnit*(cluster-1)]
            if(ord %in% c('BL','TU')) { 
                for(ii in 1:length(EVclusIncRank))
                    EVpopH[cluster==EVclusIncRank[ii] & vaccDay==Inf, vaccDay := vaccDaysLeft[ii]]
                rm(ii) }
        }

        ## ##################################################
        ## Option 2. Vacc partially vacc clusters first, then vaccinate entire unvacc clusters.
        ## vaccinate the rest of partially vaccinated clusters afterwards
        if(RCTendOption==2) {
            ## Vaccinate the rest of partially vaccinated clusters first. 
            ## First get order in which to do it, if ord!='none'
            vaccDaysLeft <- seq(min(daySeq[daySeq > endTrialDay]), 700, by = delayUnit) ## remaining vaccination days (delayUnit apart)
            ## 700 ensures that everyone eventually gets vaccinated, even if after trial stops keeping track
            if(ord == 'BL') ## reorder vaccination sequence by vaccinating highest hazard clusters as determined when trial ended
                EVclusIncRank <- EVpopH[cluster %in% vaccClusters & idByClus==1 & day == endTrialDay,
                                        rev(order(clusHaz))]
            if(ord == 'TU') {  ## reorder vaccination sequence by vaccinating highest hazard partially vacc clusters as vaccination happens
                EVclusIncRank <- NULL
                for(ii in 1:length(vaccClusters)) {
                    dd <- vaccDaysLeft[ii]
                    tmpRank <- EVpopH[cluster %in% vaccClusters & idByClus==1 & day == min(dd,maxInfectDay),
                                      rev(order(clusHaz))]
                    tmpRank <- tmpRank[!tmpRank %in% EVclusIncRank]
                    EVclusIncRank <- c(EVclusIncRank, tmpRank[1])
                }
            }
            ## Then update the vaccination days for unvaccinated individuals in partially vaccinated clusters based on these orders

            if(ord == 'none')
                EVpopH[cluster %in% vaccClusters & vaccDay==Inf, 
                       vaccDay := vaccDaysLeft[1] + delayUnit*(cluster-1)]
            if(ord == 'BL') { 
                for(ii in 1:length(EVclusIncRank))
                    EVpopH[cluster==EVclusIncRank[ii] & vaccDay==Inf, vaccDay := vaccDaysLeft[ii]]
                rm(ii) }
            if(ord == 'TU') {
                for(ii in 1:length(EVclusIncRank))
                    EVpopH[cluster==EVclusIncRank[ii] & vaccDay==Inf, vaccDay := vaccDaysLeft[ii]]
            }
            ## EVpopH[cluster %in% vaccClusters & idByClus %in% c(1,clusSize/2+1) & day == 0, ] ## check ordering works
            ## Then go back and vaccinate unvaccinated clusters

            ## print(EVpopH[idByClus %in% 1 & day == endTrialDay + (0-1)*delayUnit,  ## control group w/in each cluster
            ##        list(cluster, vaccDay, day, clusHaz, ord = order(rev(order(clusHaz))))])

            ##    EVpopH[cluster%in% vaccClusters & idByClus==1, ]

            if(length(notYetVaccClusters)>0) {
                vaccDaysLeft <- seq(min(daySeq[daySeq > EVpopH[cluster %in% vaccClusters, max(vaccDay)]]), 700, by = delayUnit)
                ## Again get order first, for totally unvaccinated clusters
                if(ord == 'BL') ## reorder vaccination sequence by vaccinating highest hazard clusters as determined when *trial ended*
                    EVclusIncRank <- EVpopH[cluster %in% notYetVaccClusters & vaccDay == Inf & idByClus==1 & day == endTrialDay,
                                            rev(order(clusHaz))]
                if(ord == 'TU') {  ## reorder vaccination sequence by vaccinating highest hazard partially vacc clusters as vaccination happens
                    EVclusIncRank <- NULL
                    for(ii in 1:length(notYetVaccClusters)) {
                        dd <- vaccDaysLeft[ii]
                        tmpRank <- EVpopH[cluster %in% notYetVaccClusters & idByClus==1 & day == dd,
                                          rev(order(clusHaz))]
                        tmpRank <- tmpRank[!tmpRank %in% EVclusIncRank]
                        EVclusIncRank <- c(EVclusIncRank, tmpRank[1])
                    }
                }
                ## Then set vaccination day for ALL within cluster (not those that
                ## haven't been already assigned, i.e. don't require vaccDay==Inf)
                if(ord == 'none')
                    EVpopH[cluster %in% notYetVaccClusters,
                           vaccDay := vaccDaysLeft[1] + delayUnit*(cluster-1)]
                if(ord %in% c('BL','TU')) { 
                    for(ii in 1:length(EVclusIncRank))
                        EVpopH[cluster==notYetVaccClusters[EVclusIncRank[ii]],
                               vaccDay := vaccDaysLeft[ii]]
                    rm(ii) }
            }
            
        }

        ## ##################################################
        ## Option 3. Continue with cluster sequence, but vaccinate full clusters. *Simultaneously* begin going
        ## back and vacc partially vacc clusters.
        if(F) {# RCTendOption==3) {

        }
    }
})

