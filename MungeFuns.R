## Construct survival data from waiting times
makeSurvDat <- function(parms,  whichDo='pop') within(parms, {
    if(verbose ==1.5) browser()
    popTmp <- get(whichDo)
    popTmp$immuneDayThink <- popTmp[,vaccDay] + immunoDelayThink ## vaccine refractory period ASSUMED in analysis
    ## pre-immunity table
    stPre <- copy(popTmp) # st = survival table
    stPre$startDay <- 0
    stPre$endDay <- stPre[, pmin(immuneDayThink, infectDay)]
    stPre$infected <- stPre[ ,as.numeric(infectDay <= immuneDayThink)]
    stPre$immuneGrp <-  0     ## immuneGrp is variable used for analysis, not omnietient knowledge of vaccination/immune status
    stPre <- stPre[,list(indiv, cluster, pair, idByClus, vaccDay, immuneDay, 
                         immuneDayThink, startDay, endDay, infectDay, infected, immuneGrp)]
    ## post-immunity table
    stPost <- copy(popTmp)[infectDay > immuneDayThink,]
    stPost$startDay <- stPost[,immuneDayThink]
    stPost$endDay   <-  stPost[,infectDay]
    stPost$infected <- 1 ## everyone gets infected eventually, but will truncate this in a separate function
    stPost$immuneGrp <- 1
    stPost <- stPost[,list(indiv, cluster, pair, idByClus, vaccDay, immuneDay, 
                           immuneDayThink, startDay, endDay, infectDay,  infected, immuneGrp)]
    nmSt <- paste0('st', sub('pop','',whichDo)) ## makes EVpop into EVst for example
    assign(nmSt, rbind(stPre, stPost)) ## combine tables
    rm(stPre, stPost, popTmp, nmSt)
}) ## careful with modifying parms, st depends on analysis a bit too (immunoDelayThink), so we can have different st for same popH

## Select subset of *s*urvival *t*able to analyze
activeFXN <- function(parms, whichDo='st') within(parms, { 
    if(verbose ==1.6) browser()
    ## for SWCT or unmatched CRCT always include all clusters in analysis because unvaccinated
    ## clusters are still considered to provide useful information from baseline
    stA <- copy(get(whichDo))
    stA$firstActive <- 0
    if(!includeAllControlPT) { ## remove person-time observed prior to post-refractory period from data
        if(trial=='CRCT' & ord!='none') ## active once anyone considered immune in matched cluster pair
            stA[, firstActive := min(immuneDayThink), by = pair]
        if(trial %in% c('RCT','FRCT')) ## active once anyone considered immune in cluster
            stA[, firstActive := min(immuneDayThink), by = cluster]
        if(trial=='SWCT') {## inactive during protective delay; active only when there exists both vacc & unvacc person-time observed
            if(remStartFin) {
                firstDayAnyoneImmune <- stA[, min(immuneDayThink)]
                lastDayAnyoneNotVacc <- stA[, max(immuneDayThink)]
                if(verbose>1.5)    {
                    print(paste0('excluding ', stA[endDay <= firstDayAnyoneImmune, sum(infectDay<Inf)], 
                                 ' infected before ', firstDayAnyoneImmune, ' days'))
                    print(paste0('excluding ', stA[startDay >= lastDayAnyoneNotVacc, sum(infectDay<Inf)], 
                                 ' infected after ', lastDayAnyoneNotVacc, ' days'))
                } 
                stA <- stA[!endDay <= firstDayAnyoneImmune] ## remove inactive observation intervals at beggining of trial
                stA <- stA[!startDay >= lastDayAnyoneNotVacc] ## remove inactive observation intervals at end of trial
                ## Left-truncate person-time before firstDay
                stA[startDay < firstDayAnyoneImmune, startDay:=firstDayAnyoneImmune]
                ## Right-truncate person-time after lastDay
                stA[endDay > lastDayAnyoneNotVacc, endDay:=lastDayAnyoneNotVacc]
                ## Now check that infectDay is still in intervals
                if(verbose>1.5) print(paste0('excluding ', 
                                             stA[, sum(infected==1 & infectDay < Inf & !(infectDay >= startDay & infectDay <= endDay))]
                                             , ' infected after ', lastDayAnyoneNotVacc, ' days'))
                stA[infected==1 & infectDay < Inf & !(infectDay >= startDay & infectDay <= endDay), infected:=0]
                rm(firstDayAnyoneImmune, lastDayAnyoneNotVacc)
            }
            if(remProtDel) { ## exclude protective delay
                ## right-truncate at vaccine date person-time intervals that starts before vaccine date and ends
                ## after vaccine date (i.e. ignore person-time within protective delay)
                stA[startDay <= vaccDay & endDay >= vaccDay, endDay := vaccDay]
                ## remove person-time completely contained within protective delay (should only remove cluster 1's 0-immunedaythink person-time, but
                ## already removed above, so redundant)
                stA <- stA[!(startDay >= vaccDay & endDay <= immuneDayThink)]
                ## Now check that infectDay is still in intervals
                if(verbose>1.5) 
                    print(paste0('excluding ', stA[, sum(infected==1 & infectDay < Inf & !(infectDay >= startDay & infectDay <= endDay))]
                                 , ' infected in protective delay'))
                stA[infected==1 & infectDay < Inf & !(infectDay >= startDay & infectDay <= endDay), infected:=0]
            }
            ## make events that happen after endDay not count
            stA[infectDay>endDay, infected:=0]
            ## run to see how person-time is distributed between cluster
            ## stA[idByClus==1, list(cluster, vaccDay, immuneDayThink, startDay, endDay)]
            ## stA[idByClus==1 & immuneGrp==1, list(cluster, vaccDay, immuneDayThink, startDay, endDay)]
            ## plotSTA(stA)
            ## stA[infectDay<Inf, list(cluster, vaccDay, immuneDayThink, startDay, endDay,immuneGrp)]
        }
    }
    stA <- stA[!endDay <= firstActive] ## remove inactive observation intervals (does nothing for SWCT)
    stA[startDay < firstActive, startDay := firstActive] ## set accumulation of person time as when the cluster/pair is active (does nothing for SWCT)
    nmStA <- paste0('stActive', sub('st','',whichDo)) ## makes EVpop into EVst for example
    assign(nmStA, stA)
    rm(stA, nmStA)
})
## p1 <- simTrial(makeParms('RCT', ord='BL', small=F), br=F)
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## s1$st[idByClus%in%1:2, list(indiv, cluster, pair, idByClus,immuneDayThink, startDay,endDay)]

plotSTA <- function(stA) {
    par(mfrow=c(2,1))
    plot(0,0, type = 'n', xlim = c(0,168), ylim = c(0,6000), bty = 'n', xlab = 'day', ylab='ID in trial', yaxt='n', main = 'infected individuals')
    stA[infected==1 & immuneGrp==0, segments(startDay, indiv, endDay, indiv, col = 'dodger blue'), cluster]
    stA[infected==1 & infectDay<Inf & immuneGrp==1, segments(startDay, indiv, min(168,endDay), indiv, col = 'red'), cluster]            
    plot(0,0, type = 'n', xlim = c(0,168), ylim = c(0,6000), bty = 'n', xlab = 'day', ylab='ID in trial', yaxt='n', main = '1st individual per cluster')
    stA[idByClus==1 & immuneGrp==0, segments(startDay, indiv, endDay, indiv, col = 'dodger blue'), cluster]
    stA[idByClus==1 & immuneGrp==1, segments(startDay, indiv, min(168,endDay), indiv, col = 'red'), cluster]
    print('# infections')
    print(stA[, sum(infected==1), immuneGrp])
    print('empirical hazard (remember declining incidence distorts this')
    print(stA[, sum(infected==1)/sum(min(endDay,168)-min(startDay,168)), immuneGrp])
}

## Restructure for GEE/GLMM with weekly observations of each cluster.
makeGEEDat <- function(parms, whichDo='popH') within(parms, {
    if(verbose ==1.7) browser()
    popHTmp <- get(whichDo)
    popHTmp$immuneDayThink <- popHTmp[,vaccDay] + immunoDelayThink ## vaccine refractory period ASSUMED in analysis
    popHTmp$infectDayRCens <- popHTmp$infectDay
    popHTmp[infectDay==Inf, infectDayRCens := NA]
    popHTmp$immuneGrp <- 0
    popHTmp[day >= immuneDayThink, immuneGrp := 1]
    popHTmp$firstActive <- 0
    if(!includeAllControlPT) { ## remove person-time observed prior to post-refractory period from data
        if(trial=='CRCT' & ord!='none') ## active once anyone considered immune in matched cluster pair
            popHTmp[, firstActive := min(immuneDayThink), by = pair]
        if(trial %in% c('RCT','FRCT')) ## active once anyone considered immune in cluster
            popHTmp[, firstActive := min(immuneDayThink), cluster]
    }
    popHTmp$active <- popHTmp[,day>=firstActive]
    clusD <- popHTmp[, list(cases = sum(!is.na(infectDayRCens)), atRisk = length(idByClus)), 
                     list(cluster, day, active, immuneGrp, vaccDay, immuneDayThink)]
    clusD <- clusD[order(cluster, day)]
    clusD[, atRisk := atRisk - c(0L, cumsum(cases[-length(cases)])), cluster]
    clusD <- clusD[active==T]
    if(!includeAllControlPT & trial == 'SWCT') {
        if(remStartFin) {
            firstDayAnyoneImmune <- popHTmp[, min(immuneDayThink)]
            lastDayAnyoneNotVacc <- popHTmp[, max(immuneDayThink)] - 1
            clusD <- clusD[!day < firstDayAnyoneImmune] ## remove inactive observation intervals at beggining of trial
            clusD <- clusD[!day > lastDayAnyoneNotVacc] ## remove inactive observation intervals at end of trial
            rm(firstDayAnyoneImmune, lastDayAnyoneNotVacc)
        }
        if(remProtDel) { ## exclude protective delay
            clusD <- clusD[!(day >= vaccDay & day < immuneDayThink)]
        }
    }
    ## plotClusD(clusD) ## to check result
    nmSt <- paste0('clusDat', sub('popH','',whichDo)) ## makes EVpop into EVst for example
    assign(nmSt, clusD) ## combine tables
    rm(popHTmp, clusD, nmSt)
})

plotClusD <- function(clusD) {
    plot(0,0, type = 'n', xlim = c(0,168), ylim = c(0,20), bty = 'n', xlab = 'day', ylab='ID in trial', yaxt='n', main = 'infected individuals')
    clusD[immuneGrp==0, segments(day, cluster, day+7, cluster, col = 'dodger blue'), cluster]
    clusD[immuneGrp==1, segments(day, cluster, day+7, cluster, col = 'red'), cluster]    
    print('# infections')
    print(clusD[, sum(cases), immuneGrp])
    print('empirical hazard (remember declining incidence distorts this')
    print(clusD[, sum(cases)/(sum(atRisk)*7), immuneGrp])
}

## Take a survival data from above function and censor it by a specified time in months
censSurvDat <- function(parms, censorDay = parms$maxDurationDay, whichDo = 'stActive') with(parms, {
    if(verbose==2.7) browser()
    stTmp <- copy(get(whichDo))
    intervalNotStarted <- stTmp[,startDay] > censorDay
    stTmp <- stTmp[!intervalNotStarted,] 
    noInfectionBeforeCensor <- stTmp[,endDay] > censorDay
    stTmp[noInfectionBeforeCensor, infected:=0]
    stTmp[noInfectionBeforeCensor, endDay:=censorDay]
    stTmp[,perstime := (endDay-startDay)]
    stTmp <- stTmp[perstime > 0,] 
    return(stTmp)
})

summTrial <- function(st) list(summarise(group_by(st, cluster), sum(infected))
                               , summarise(group_by(st, cluster, immuneGrp), sum(infected))
                               , summarise(group_by(st, immuneGrp), sum(infected))
                               )
