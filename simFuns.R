require(blme); require(survival); require(coxme); require(data.table); require(parallel); require(dplyr); require(msm)
load('data/createHT.Rdata')

yearToDays <- 1/365.25
monthToDays <- 1/30
trialTypes <- c('RCT','FRCT','SWCT','CRCT')
makeParms <- function(
    trial='RCT'
  , trialStartDate = '2015-02-01' ## converted to date below    
  , numClus=20, clusSize=300
  , delayUnit = 7 ## logistically imposed interval in between each new cluster receiving vaccination
  , ord = 'none' ## order clusters' receipt of vaccination ('none', by baseline visit 'BL', by time-updated 'TU' last interval incidence)
  , hazType =  'SL' ## use hazards from "SL" or "Phenom"enologically driven hazards
  , HazTrajSeed = NA ## allows keeping all forecasted hazard trajectories identical between simulations (only demograhic stochasticity)
  , nbsize = .8 ## for above
  , propInTrial = .03 ## for above
  , mu=.03 * yearToDays ## mean hazard in all participants at baseline
  , cvClus=1 ##  variance in cluster-level hazards for gamma distribution
  , cvClusTime=1 ##  temporal fluctuation variance in cluster-level hazards around smooth trajectories for gamma distribution 
  , sdLogIndiv = 1 ## variance of lognormal distribution of individual RR within a hazard (constant over time, i.e. due to job)
  , vaccEff = .8
  , pSAE = 10^-4
  , maxDurationDay = 7*24 ## maximum duration end of trial (24 weeks default; 6 months) (trial can stop early though; e.g. endTrialDay)
  , trackUntilDay = 2*maxDurationDay ## how long to track infections for calculating incidence averted & equipoise
    ## Sequential Design info
  , gs=FALSE
  , gsDesArgs = list(k=3, test.type=2, alpha=0.025, beta=0.1, timing = seq(0,1, l = 4)[-1]) ## total number of analyses
  , infoCalcArgs = list(assumedVaccEff = .7, assumedCumHazard = 10^-4, beta = .1) ## number of events to power trial asymptotes at low cumulative hazards (e.g. 36 for vaccEff=.7)
  , doCFs=F ## do coutnerfactual simulations
  , numCFs = 1
    ## , seqType='MaxDur' ## MaxDur or MaxInfo, i.em. spend alpha based on # events up until endtime, then spend rest of alpha
  , immunoDelay = 21 ## delay from vaccination to immunity
  , immunoDelayThink = immunoDelay ## delay from vaccination to immunity used in analysis (realistically would be unknown)
  , weeklyDecay=.9, cvWeeklyDecay=1 ## log-normally distributed incidence decay rates (set var = 0 for constant)
  , hazIntUnit = 7 ## interval between discrete changes in hazard
  , reordLag = 14 ## how long ago's hazard to use when deciding this week's time-updated vaccination sequence
  , includeAllControlPT = F ## include person-time from controlled trials before end of
    ## vaccination refractory period? in SWCT also only includes
    ## person-time when trial has both protected & unprotected individuals
  , remProtDel = T## remove protective delay PT from SWCT
  , remStartFin = F ## remove start
  , RCTendOption = 2        ## order to vaccinate unvaccinated invididuals when an RCT ends, see EndTrialFuns.R
  , instVaccDelay = 7 ## delay til instant vacc of everyone after trial ends in trials where delayUnit=0 otherwise
  , small=F ## do a small trial for illustration
  , StatsFxns = c('doCoxME')##,'doGLMFclus','doGMMclus','doGLMclus','doRelabel','doBoot')
  , nboot = 200 ## bootstrap samples
  , doInf = T
  , verbose = 0
  , dontReordForPlot = F
){
    if(small) {
        numClus <- 4
        clusSize <- 4
    }
    trialStartDate <- as.Date(trialStartDate)
    if(maxDurationDay < delayUnit*numClus) stop('maxDurationDay too short. Need enough time to rollout vaccines to all clusters')
    if(gsDesArgs$k>1 & length(StatsFxns)>1) stop('Cannot apply multiple statistical functions for simulations of sequential designs since trial progression is conditional on analysis')
    if(trial=='FRCT') delayUnit <- delayUnit/2 ## rolling out vaccines as quickly as you would if you were vaccinating whole clusters
    if(gs) {
        gsBounds <- do.call(gsDesign, gsDesArgs)
        maxInfo <- with(infoCalcArgs, { ## total number of events expected need to power the trial with 1-beta power for an assumedVaccEff 
            cumHazs <- assumedCumHazard*c(1, (1-assumedVaccEff))
            sampSizes <- nBinomial(p1 = cumHazs[1], p2 = cumHazs[2], outtype=2, beta = beta)
            return(sum(sampSizes*cumHazs))
        })
    }
    ## if(!is.null(vaccProp)) {
    ##     vaccEff <- vaccProp[numInBatch, vaccEff]
    ##     pSAE <- vaccProp[numInBatch, pSAE]
    ## }
    return(as.list(environment()))
}

## Make a trial population with a given number of clusters of a given size. Put the people in
## clusters, give them individual IDs and also id # within cluster
makePop <- function(parms=makeParms()) within(parms, {
    pop <- data.table(indiv=as.factor(1:(numClus*clusSize))
                    , cluster=as.numeric(gl(n=numClus, k=clusSize))
                    , idByClus = rep(1:clusSize, numClus)
                      )
})

## reparameterize a gamma by mean/var to simulate spatial variation in underlying hazards (change
## later to something more reasonable or based on real data)
reParmRgamma <- function(nn, mean, cv) {
    if(cv > 0) {
        var <- (cv*mean)^2 ## cv = sd/mean, so sd = cv*mean, so var = (cv*mean)^2
        theta <- var/mean
        k <- mean/theta
        rgamma(nn, shape = k, scale = theta)
    }else{ ## no variance
        mean
    }
}

## Phenomenological Hazard Trajectories to explore exact cause of elevated Type I errors
createHazTraj_Phenom <- function(parms) within(parms, {
    if(verbose>10) browser()
    baseClusHaz <- reParmRgamma(numClus, mean = mu, cv = cvClus) ## gamma distributed baseline hazards
    wd <- weeklyDecay
    inverted <- F
    if(weeklyDecay>1) { ## so we can use logit transform for growing epidemics too
        wd <- 1/weeklyDecay
        inverted <- T
    }
    dailyDecayRates <- inv.logit(rnorm(numClus, mean = logit(wd^(1/7)), sd = cvWeeklyDecay*logit(wd)))
    if(inverted) dailyDecayRates <- 1/dailyDecayRates
    hazT <- data.table(day = rep(daySeq, each = numClus), cluster = rep(1:numClus, length(daySeq)), clusHaz = 0)
    cHind <- which(names(hazT)=='clusHaz')
    ## mean cluster hazard trajectory
    for(ii in 1:numClus) hazT[which(hazT[,cluster]==ii), clusHaz := baseClusHaz[ii]*dailyDecayRates[ii]^day]
    ## negative binomial variation around smooth declines
    hazT[, clusHaz := reParmRgamma(nrow(hazT), mean = clusHaz, cv = cvClusTime)]
    rm(ii, cHind, baseClusHaz)
})

## Hazard Trajectories from SL district-level incidence (calls fits peviously made)
createHazTraj_SL <- function(parms) within(parms, {
    if(verbose>10) browser()
    hazT <- data.table(createHazTrajFromSLProjection(fits, trialStartDate = trialStartDate, HazTrajSeed=HazTrajSeed,
                                                     nbsize = nbsize, propInTrial = propInTrial, verbose=verbose,
                                                     clusSize = clusSize, numClus = numClus, weeks = T))
    hazT <- hazT[day %in% daySeq]
    hazT[clusHaz==0, clusHaz := 10^-8] ## to stablize things
    setcolorder(hazT, c('day','cluster','clusHaz','Date'))
})

setClusHaz <- function(parms) {
    parms <- within(parms, { ## add daySeq
        daySeq <- seq(-hazIntUnit*ceiling(reordLag/hazIntUnit),trackUntilDay,by=hazIntUnit) })
    tempCHTfxn <- get(paste0('createHazTraj_', parms$hazType)) ## get chosen function
    return(tempCHTfxn(parms))
}
## parms <- makePop()
## setClusHaz(within(parms, {hazType='SL'}))$hazT
## setClusHaz(within(parms, {hazType='Phenom'}))$hazT


## Set cluster- and individual-level hazards, with cluster means changing over time and individual
## RR around cluster mean constant
setIndHaz <- function(parms=makePop()) within(parms, {
    if(verbose>10) browser()
    ## give every individual a lognormally distributed relative risk
    pop$indivRR <- rlnorm(numClus*clusSize, meanlog = 0, sdlog = sdLogIndiv)
    ## create popH which has weekly hazards for all individuals
    popH <- arrange(merge(pop, hazT, by='cluster', allow.cartesian=T),day)
    popH$indivHaz <- popH[, clusHaz*indivRR]
    daySeq <- daySeq[daySeq>=0] ## don't do anything before zero, just stored hazard in popHearly for ordering
    daySeqLong <- seq(0,trackUntilDay+1000,by=hazIntUnit) ## to avoid problems later
    popHearly <- copy(popH)
    popH <- popH[day >= 0]
    rm(hazT)
})
## setHazs(makePop(makeParms(weeklyDecay=1, weeklyDecayVar=0)))$popH[cluster==1,]
## setHazs(makePop(makeParms(weeklyDecay=.9, weeklyDecayVar=0)))$popH[cluster==1,]

reordPop <- function(parms) { ## wrapper around other functions below
    reordFXN <- get(paste0('reord',parms$trial))
    parms <- reordFXN(parms)
    within(parms, { ## return parms
        if(verbose>10) browser()
        popH[, cluster:=clusIncRank[popH[, cluster]]]
        popH <- arrange(popH, day, cluster)
        ## reset indices so they're ordered again by vaccination sequence
        popH[,indiv:= rep(1:(numClus*clusSize), popH[, length(unique(day))])] 
        rm(clusIncRank)
        popH$pair <-  NA ## pairs matched for randomization (if matching)
        if(trial=='CRCT' & ord!='none') { ## only paired clusters exist in matched CRCT
            popH$pair <- popH[, cluster %% (numClus/2)]
            popH[pair==0, pair:=numClus/2]
        }
    })
}

reordSWCT <- reordFRCT <- reordRCT <- function(parms) within(parms, {
    if(verbose>10) browser()
    if(ord=='none') { ## should already be random but do it agan for good measure (debugging randomziation)
        clusIncRank <- sample(1:numClus, numClus, replace = F)
        if(dontReordForPlot) clusIncRank <- 1:20
    }
    if(ord=='BL') { ## if vaccinating highest incidence clusters first (i.e. *NOT* SW randomization of vaccination sequence)
        clusIncRank <- popH[idByClus==1 & day == 0,order(rev(order(clusHaz)))]
    }
    if(ord=='TU') { ## time-updated sequencing, each week vaccinate highest incidence cluster from reordLag days ago
        clusIncRank <- NULL
        for(ii in 1:numClus) { ## for each vaccination day
            dd <- daySeq[ii]
            updatingOrder <- 1:numClus > ii-1 ## i.e. on 3rd day of vaccination, only updating the 3rd vaccination sequence
            currentRank <- popHearly[idByClus==1 & day == dd - reordLag, rev(order(clusHaz))] ## current cluster hazard ordering
            currentRank <- currentRank[!currentRank %in% clusIncRank] ## remove clusters already vaccinated
            clusIncRank <- c(clusIncRank, currentRank[1])
        }
        clusIncRank <- order(clusIncRank)
        rm(currentRank,updatingOrder,dd,ii)
    }
})

reordCRCT <- function(parms) within(parms, {
    if(verbose>10) browser()
    if(ord=='none') { ## should already be random but do it agan for good measure (debugging randomziation)
        clusIncRank <- sample(1:numClus, numClus, replace = F)
    }
    if(ord=='BL') { ## if matching clusters on incidence, then randomizing within pairs, then vaccinating highest incidence first
        if(numClus %% 2 == 1) stop("need even # of clusters")
        ## order clusters by mean hazard, in pairs (each row)
        clusIncRank <- matrix(rev(order(popH[day==0 & idByClus==1, clusHaz])), nc = 2, byrow=T) 
        whichVacc <- sample(1:2, nrow(clusIncRank), replace=T) ## which of each pair to vaccinate
        for(ii in 1:nrow(clusIncRank)) if(whichVacc[ii]==2) clusIncRank[ii,1:2] <- clusIncRank[ii,2:1] ## make first column vaccinated group
        clusIncRank <- c(clusIncRank[,1], clusIncRank[,2]) ## reorder of clusters
        clusIncRank <- order(clusIncRank)
        rm(whichVacc, ii)
    }
    ## time-updated sequencing, each week select the two highest incidence pairs that have yet to be
    ## randomized, and randomize one of them to vaccination
    if(ord=='TU') { 
        if(numClus %% 2 == 1) stop("need even # of clusters")
        numPairs <- numClus/2
        clusIncRank <- NULL
        for(ii in 1:numPairs) { ## for each vaccination day
            dd <- daySeq[ii]
            notRandomized <- (1:numClus)[! 1:numClus %in% as.vector(clusIncRank)] ## haven't already been randomized
            currIncOrd <- notRandomized[rev(order(popHearly[day==dd - reordLag & idByClus==1 & cluster %in% notRandomized, clusHaz]))]
            currentRank <- matrix(currIncOrd, nc = 2, byrow=T)[1,] ## pick top row of paired matrix
            if(rbinom(1,1,.5)) currentRank <- rev(currentRank) ## deterine which of each pair to vaccinate (1st column)
            clusIncRank <- rbind(clusIncRank, currentRank) ## add pair
        }
        clusIncRank <- as.numeric(c(clusIncRank[,1], clusIncRank[,2])) ## vaccinated, control
        clusIncRank <- order(clusIncRank) ## get reording for next line
        rm(notRandomized,dd,ii,numPairs,currIncOrd,currentRank)
    }
})

## p1 <- setHazs(makePop(makeParms('CRCT', clusSize=2, weeklyDecay=.9, weeklyDecayVar=.3, ord='BL')))
## p1 <- reordPop(p1)
## p1$popH
## p1$popH[idByClus==1 & day==0 & cluster <=numClus/2, order(clusHaz)]

setVaccDays <- function(parms) { ## wrapper around other functions below
    setVaccFXN <- get(paste0('set',parms$trial,'vaccDays'))
    parms <- setVaccFXN(parms)
    return(parms)
}
setSWCTvaccDays <- function(parms, whichDo='pop') within(parms, {
    tmpstrg <- paste0(whichDo, 'H')
    tmpH <- copy(get(tmpstrg))
    tmpH$vaccDay <- delayUnit*(tmpH[,cluster]-1)
    assign(tmpstrg, tmpH)
    rm(tmpH, tmpstrg)
})
setFRCTvaccDays <- setRCTvaccDays <- function(parms) within(parms, { ## assuming same speed rollout as SWCT (unless FRCT)
    popH$vaccDay <- Inf ## unvaccinated
    popH[idByClus > clusSize/2 , vaccDay := delayUnit*(cluster-1)] ## half get vaccinated in each cluster, 1 per interval
})
setCRCTvaccDays <- function(parms) within(parms, {
    popH$vaccDay <- Inf
    popH[cluster <= numClus/2, vaccDay := delayUnit*(cluster-1)] ## first half of clusters (1 per pair) get vaccinated in sequence
})
## To check ordering works
## p1 <- setHazs(makePop(makeParms(clusSize=2, weeklyDecay=.9, weeklyDecayVar=.2, ord='TU', trial='SWCT',small=T)))
## setVaccDays(p1)$popH[,list(cluster,clusHaz, day,vacc,immune)]
## p1 <- setVaccDays(p1)
## p1$popH[idByClus==1,list(cluster,clusHaz, day,vacc,immune)]
setImmuneDays <- function(parms, whichDo='pop') within(parms, {
    tmpstrgH <- paste0(whichDo, 'H')
    tmpH <- copy(get(tmpstrgH))
    tmpH$immuneDay <- tmpH[,vaccDay] + immunoDelay ## vaccine refractory period
    tmpH$vacc <- tmpH[, day>=vaccDay]
    tmpH$immune <- tmpH[, day>=immuneDay]
    ## reset pop to refrence data table after reordering and then assignment of vaccday stuff
    tmp <- select(tmpH[day==0], indiv, cluster, pair, idByClus, indivRR, vaccDay, immuneDay) 
    assign(tmpstrgH, tmpH)
    assign(whichDo, tmp)
    rm(tmpH, tmpstrgH)
})

## Simulate infections. Takes popH for initial simulation, or EVpopH for end trial vaccination version (requires startInf)
simInfection <- function(parms, whichDo='pop', startInfectingDay = 0, cfNum=1, RNGseed = NULL) ## startInf can be set to endTrialDay
    within(parms, { 
        if(verbose>10 | verbose == 9.89) browser()
        tmp <- get(whichDo)
        tmpH <- get(paste0(whichDo,'H'))
        if(startInfectingDay==0) tmp$infectDay <- tmpH$infectDay <- Inf ## otherwise it's already got some infection data in it
        tmpH[infectDay > startInfectingDay, infectDay := Inf] ## redoing post endDay stuff with additional folks vacc
        tmp[infectDay > startInfectingDay, infectDay := Inf] ## redoing post endDay stuff with additional folks vacc
        alreadyInfected <- NULL
        ## RNG seed control
        ## if(whichDo %in% c('pop','NTpop','VRpop')) 
        ##     assign(paste0('simInfSeed',whichDo), .GlobalEnv$.Random.seed) ## saved for use elsewhere later
        if(!is.null(RNGseed))
            .GlobalEnv$.Random.seed <- RNGseed 
        ## Loop thru infections: infection day is beginning of each hazard interval + exponential waiting time
        for(dd in daySeq[daySeq>=startInfectingDay]) { 
            alreadyInfected <- tmpH[infectDay!=Inf, indiv] ## don't reinfect those already infected
            tmpH[day==dd & !indiv %in% alreadyInfected, 
                 infectDay := dd + rexp(length(indiv), rate = indivHaz*ifelse(immune, 1-vaccEff, 1))] 
            tmpH[day==dd & !indiv %in% alreadyInfected & infectDay > dd + hazIntUnit, 
                 infectDay := Inf] ## reset if it goes into next hazard interval
        }
        ## copy infection days to pop, to use in analysis
        indivInfDays <- tmpH[infectDay!=Inf & infectDay > startInfectingDay, list(indiv,infectDay, vaccDay, indivRR)]
        indivInfDays <- arrange(indivInfDays, indiv)
        tmp[indiv %in% indivInfDays[,indiv], infectDay:= indivInfDays[,infectDay]]
        ## SAEs
        tmp[, SAE:= as.integer(rbinom(length(indiv), 1, pSAE))]
        tmp[vaccDay==Inf, SAE:=0] ## can't have SAE if was not vaccinated
        tmp[infectDay<vaccDay, SAE:=0] ## can't have SAE if got Ebola before vaccination (excluded from vaccination in any case)
        tmpH$SAE <- 0
        setkey(tmpH, indiv, day) ## fastest way to index by multiple columns
        tmpH[tmp[SAE==1, list(indiv, vaccDay)], SAE:=1]
        if(whichDo %in% c('NTpop','VRpop')) {## output summary of the population's infections & SAEs
            indivEventDays <- tmp[infectDay!=Inf | SAE==1, list(indiv,infectDay, vaccDay, SAE, indivRR)]
            if(nrow(indivEventDays)==0) ## still need a row in this data.table, otherwise messes up summaries later in finInfo
                indivEventDays <- data.table(indiv=NA, infectDay=Inf, vaccDay=Inf, SAE=0, indivRR=NA)
        }
        ## Assign back to global variables
        assign(whichDo, tmp)
        assign(paste0(whichDo,'H'), tmpH)
        rm(tmp, tmpH, dd,alreadyInfected, indivInfDays)
    })

## p1 <- setHazs(makePop(makeParms(clusSize=300, numClus=20, weeklyDecay=.9, weeklyDecayVar=0, ord='BL')))
## p1 <- reordPop(p1)
## p1 <- simInfection(p1)
## head(p1$pop[infectDay!=Inf, list(cluster, immuneDay, infectDay)],100)

## simulate whole trial and return with all parameters used
simTrial <- function(parms=makeParms(), seed = NULL) {
    if(!is.null(seed)) set.seed(seed) ## for plotting comparisons between hazard trends w/ & w/o fluctuationrs
    parms <- makePop(parms) ## make population
    parms <- setClusHaz(parms) ## set cluster-level hazards
    parms <- setIndHaz(parms) ## set individual-level hazards
    parms <- reordPop(parms) ## reorder vaccination sequence by incidence (if applicable), doReord off for plotting only
    parms <- setVaccDays(parms) ## set vaccination days
    if(parms$doInf) {
        parms <- setImmuneDays(parms) ## set vaccination days
        parms <- simInfection(parms) ## simulate infection
    }
    return(parms)
}
## simTrial(b=F)
## simTrial(makeParms(small=T),F)

subsArgs <- function(parms, fxn) parms[names(parms) %in% names(formals(fxn))] ## get parameters necessary for a fxn

## For setting up batch runs (see *MK.R files)
addParm <- function(x, parmsMat,ii, index='rcmdbatch') {
    for(pp in 1:length(parmsMat)) {
        tempP <- as.data.frame(parmsMat)[,pp]
        isch <- !is.numeric(tempP[1])
        parmAdd <- tempP[parmsMat[,get(index)]==ii]
        addStrg <- paste0(" ", names(parmsMat)[pp], "=", "\""[isch], parmAdd, "\""[isch])
        x <- paste0(x, addStrg)
    }
    return(x)
}


