
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
        EVpop <- copy(pop) ## copy other infection data in to use for pre-end trial time
        EVpop$vaccDay <- EVpopH[day==0, vaccDay]
        EVpop$vaccDay <- EVpopH[day==0, immuneDay]
        EVpop$vaccDay <- EVpopH[day==0, infectDay]
    })
    ## do infection process again post-end of trial
    if(with(parms, vaccEffEst['p'] < .05 & vaccEffEst['lci']>0)) { 
        parms <- simInfection(parms, whichDo = 'EVpop', startInfectingDay = parms$endTrialDay)
    }
    ## calculate # infected for various permutations
    browser()
    parms <- makeSurvDat(parms, whichDo='EVpop')
    parms <- activeFXN(parms, whichDo = 'EVst')
    ## HERE**** haven't gotten the # cases calculated well yet
    summTrial(censSurvDat(parms, parms$maxInfectDay, 'stActive'))[[3]]
    summTrial(censSurvDat(parms, parms$maxInfectDay, 'st'))[[3]]
    summTrial(censSurvDat(parms, parms$maxInfectDay, 'EVst'))[[3]]
    summTrial(censSurvDat(parms, parms$maxInfectDay, 'EVstActive'))[[3]]


    return(parms)
}

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
