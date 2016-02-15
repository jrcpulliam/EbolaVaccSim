popNms <- c("Spop", "SpopH", "SpopEvents", "SEVpopEvents")
dparms0 <- c('trial','gs','doSL','propInTrial','nbsize', 'avHaz',
                            'ord','reordLag','delayUnit' ,'immunoDelay','trialStartDate', 'HazTrajSeed'                          
                                   )

quantcut <- function(x, qs = seq(0,1, l = 6)) {
    brks <- unique(quantile(x, qs))
    if(length(brks)==1) brks <- brks*c(.99,1/.99) ## in case it's one number, just make an interval that includes it
    return(as.numeric(cut(x, brks, include.lowest = T)))
}


extractOneSim <- function(fileNm ## prepare each simulation for binding into a large data.table
                        , dparms = dparms0
                        , verbose = 0
                        , verbFreq = 100
                        , indivLev = F 
                          ) {
    load(fileNm)
    riskStratList <-  finInfoList <- finModList <- parmsList <- list(NULL)
    if(exists('sim')) {
        if(verbose>2) browser()
        nbatch <- as.numeric(gsub("[^0-9]", "", fileNm))
        if(nbatch %% verbFreq == 0) print(nbatch)
        sim$parms[['trialStartDate']] <- as.character(sim$parms[['trialStartDate']])
        parmsList <- data.table(nbatch = nbatch, t(unlist(sim$parms[dparms])))
        finModList <- data.table(nbatch = nbatch, sim$sim$finMods)
        finInfoList <- data.table(nbatch = nbatch, sim$sim$finInfo)
        essence <- list(parms=parmsList, finTrials=finModList, finInfo=finInfoList)
        if(indivLev) {
            riskStratList <- sim$sim$Spops
            riskStratList <- within(riskStratList, {
                setkey(Spop, simNum, indiv, cluster)
                setkey(SpopH, simNum, cluster)
                setkey(SpopEvents, simNum, indiv)
                notCF <- nrow(SEVpopEvents)>0
                Spop <- merge(Spop, SpopEvents, by = c('sim','simNum','indiv'))
                if(notCF) {
                    Spop <- merge(Spop, SEVpopEvents, suffixes = c('','_EV'), by = c('sim','simNum','indiv'))
                }else{ ## or make empty columns
                    Spop <- within(Spop, { vaccDay_EV <- infectDay_EV <- SAE_EV <- NA})
                }
                Spop <- data.table(nbatch = nbatch, Spop)
                SpopH <- data.table(nbatch = nbatch, SpopH) 
                rm(notCF)
            })
            riskStratList$SpopEvents <- riskStratList$SEVpopEvents <- NULL
            essence <- c(essence, riskStratList)
        }
        return(essence)
    }
}

## Extract from simulations from multiple cores
extractSims <- function(thing
                      , dparms = dparms0
                      , doSave=T
                      , verbose = 0
                      , maxbatches = NA
                      , nbatchDo=NA
                      , mc.cores = 48
                      , indivLev = F 
                        ) {
    if(verbose==1) browser()
    batchdirnm <- file.path('BigResults',thing)
    fls <- list.files(batchdirnm, pattern=thing, full.names = T)
    nbatch <- length(fls)
    if(!is.na(maxbatches)) {
        nbatch <- maxbatches
        fls <- fls[1:nbatch]
    }
    if(!is.na(nbatchDo[1])) {
        ## batchdirnm <- file.path('BigResults',thing)
        ## flstmp <- list.files(batchdirnm, pattern=thing)
        fns <- as.numeric(gsub("[^0-9]", "", fls)) ##as.numeric(sub('.Rdata','', sub(thing,'',flstmp)))
        fls <- fls[fns %in% nbtd]
    }
    print(paste0('extracting from ', length(fls), ' files'))
    resList <- list()
    tmp <- mclapply(fls, extractOneSim, indivLev = indivLev, verbose=0, mc.cores=mc.cores)
    ## tmp <- mclapply(fls[[5340]], extractOneSim, indivLev = indivLev, verbose=1, mc.cores=mc.cores)
    length(resList) <- length(tmp[[1]])
    names(resList) <- names(tmp[[1]])
    for(vv in names(resList)) {
        resList[[vv]] <- rbindlist(lapply(tmp, function(x) {x[[vv]]}), fill=T)
        if(vv!='parms') setkey(resList[[vv]], nbatch ,simNum)
    }
    if(doSave) {
        save(resList, file=file.path('BigResults', paste0(thing, '.Rdata')))
    }
    return(resList)
}

procResList <- function(resList, verbose = 0, doSave=T) {
    resList <- within(resList, {
        if(verbose>0) browser()
        ## Coverage
        finTrials[, cvr := lci < vaccEff & uci > vaccEff]
        ## Bias, must be done on RH/(RH+1) scale to deal with Inf & 0 RH's
        finTrials$RH <- finTrials[, 1-mean]
        finTrials$PHU <- finTrials[, RH/(RH+1)] ## proportion of hazard unavoidable even with vaccination
        finTrials[RH==Inf, PHU:=1] ## otheriwse gives NaN for Inf/(Inf+1)
        finTrials[, stopped:=vaccGood|vaccBad]
        ## finTrials[,list(vaccEff,mean,PHU)] ## NEED TO AVERAGE BIAS ON PHU scale 
        print("Distribution of finTrials$err") ## -1 couldn't fit coxme, so fit coxph instead, 1 couldn't fit anything at all
        print(finTrials[, table(err)])         ## >1 is # of times got NA for permuted/bootstrapped statistic

        ## merge all trial-level info
        finit <- merge(parms, merge(finInfo, finTrials, all.x=T), all.y=T, by = 'nbatch')

        finit[order(gs), list(tcalMean = mean(tcal), power = mean(vaccGood[vaccEff>0]), length(tcal[vaccEff>0])),
              list(trial, gs, ord, delayUnit, propInTrial, trialStartDate, avHaz, cat)]
        cols <- c('trial', 'gs', 'ord', 'delayUnit', 'propInTrial', 'trialStartDate', 'avHaz', 'cat')
        finit[, (cols):=lapply(.SD, as.factor), .SDcols=cols]
        rm(cols)
        ## make sure all simulations completed (2040 of each)
        ## nsms <- finit[, .N, list(trial, gs, ord, delayUnit, propInTrial, trialStartDate, avHaz, cat)]
        ## cols <- c('trial', 'gs', 'ord', 'delayUnit', 'propInTrial', 'trialStartDate', 'avHaz', 'cat')
        ## nsms[, (cols):=lapply(.SD, as.factor), .SDcols=cols]
        ## summary(nsms[N<2040]) ## SWCT's are missing results, working on this

        finit[, posv:= vaccEff>0] ## positive vaccine efficacy simulations
        ## names(finit)
        finit[,list(trial, simNum, cat, caseTot)]

        ## Figure out how to match trials
        finit[,.N,list(nbatch, simNum, cat)] ## full unique
        finit[trial=='NT',.N,list(propInTrial, trialStartDate, avHaz, simNum)] ## unique CF
        finit[trial=='VR',.N,list(propInTrial, trialStartDate, avHaz, simNum)] ## unique CF
        finit[cat=='allFinalEV',.N,list(propInTrial, trialStartDate, avHaz, simNum)] ## matched scenarios (i.e. 7 factuals/2cfs)

        finit[, c('caseTotNT','caseTotVR'):=as.numeric(NA)]
        ## EV
        finit[cat=='allFinalEV', caseTotNT:=caseTot[trial=='NT'], list(propInTrial, trialStartDate, avHaz, simNum)]
        finit[cat=='allFinalEV', caseTotVR:=caseTot[trial=='VR'], list(propInTrial, trialStartDate, avHaz, simNum)]
        ## noEV
        finit[(trial %in% c('RCT','SWCT') & cat=='allFinal_noEV')|(trial=='NT' & cat=='allFinalEV'), caseTotNT:=caseTot[trial=='NT'],
              list(propInTrial, trialStartDate, avHaz, simNum)] #     
        finit[(trial %in% c('RCT','SWCT') & cat=='allFinal_noEV')|(trial=='VR' & cat=='allFinalEV'), caseTotVR:=caseTot[trial=='VR'],
              list(propInTrial, trialStartDate, avHaz, simNum)] #     
        ## inf Avert
        finit[, infAvert := caseTotNT - caseTot]
        finit[, list(propInTrial, trialStartDate, avHaz, simNum, trial, gs, ord, cat, caseTot, caseTotNT, caseTotVR, infAvert)] #     
        finit[, infAvertProp := infAvert/caseTotNT]
        finit[, infAvertableProp := infAvert/(caseTotNT-caseTotVR)]
        finit[, lab:=trial]
        finit[trial=='RCT' & gs==T, lab:=paste0(lab,'-gs')]
        finit[trial=='RCT' & ord=='TU', lab:=paste0(lab,'-rp')]
        finit[,lab:=as.factor(lab)]
    })
    if(doSave) {
        save(resList, file=file.path('BigResults', paste0(thing, '.Rdata')))
    }
    return(resList)
}

## Make a data table of risk-strata level info on infections spent/averted & power contributed by
## that group within the trial
makeInfPow <- function(resList, doSave=T, verbose = 0) within(resList, {
    if(verbose>0) browser()

    ## unique identifier for trial population & for parameter set
                   
    
    ##     if(doSave) {
    ##         save(resList, file=file.path('BigResults', paste0(thing, '.Rdata')))
    ##     }

})



##     tmp <- copy(get(paste0('strat_',whichDo)))
##     setnames(tmp, whichDo, 'hazCat')
##     if(!is.na(nbatchDo[1])) tmp <- tmp[nbatch %in% nbatchDo] ## do subset
##     ## make sure all combinations are filled in by doing a cross-join for all categories of batch, simnum, hazard cat & vaccdayrand
##     shc <- CJ.dt(unique(tmp[,list(nbatch, simNum)]), CJ(hazCat = tmp[,unique(hazCat)], vaccDay_noEV = tmp[,unique(vaccDay_noEV)]))
##     setnames(shc, 'hazCat', whichDo)
##     setnames(tmp, 'hazCat', whichDo)

##     setkeyv(shc, c('nbatch', 'simNum', whichDo, 'vaccDay_noEV'))
##     setkeyv(tmp, c('nbatch', 'simNum', whichDo, 'vaccDay_noEV'))
##     shc <- merge(shc, tmp, all=T) ##
##     rm(tmp); gc()

##     shc[is.na(Ninf_noEV), Ninf_noEV:=0] ## replace NA's with 0s
##     shc[is.na(Nsae_noEV), Nsae_noEV:=0] ## replace NA's with 0s
##     shc[is.na(N), N:=0]
##     ## replace NAs w 0s only if a factual (otherwise should be NA
##     shc[!nbatch %in% parms[trial %in% c('NT','VR'), nbatch] & is.na(Ninf_EV), Ninf_EV:=0] 
##     shc[!nbatch %in% parms[trial %in% c('NT','VR'), nbatch] & is.na(Ninf_noEV), Ninf_noEV:=0] 
##     shc[!nbatch %in% parms[trial %in% c('NT','VR'), nbatch] & is.na(Nsae_EV), Nsae_EV:=0]
##     shc[!nbatch %in% parms[trial %in% c('NT','VR'), nbatch] & is.na(Nsae_noEV), Nsae_noEV:=0]
##     shc[nbatch %in% parms[trial %in% c('NT','VR'), nbatch] & is.na(Ninf), Ninf:=0]
##     shc[nbatch %in% parms[trial %in% c('NT','VR'), nbatch] & is.na(Nsae), Nsae:=0]
##     shc[nbatch %in% parms[trial %in% c('NT','VR'), nbatch], Ninf_EV:=Ninf]
##     shc[nbatch %in% parms[trial %in% c('NT','VR'), nbatch], Ninf_noEV:=Ninf]
##     shc[nbatch %in% parms[trial %in% c('NT','VR'), nbatch], Nsae_EV:=Nsae]
##     shc[nbatch %in% parms[trial %in% c('NT','VR'), nbatch], Nsae_noEV:=Nsae]
##     shc$Ninf <- NULL
##     shc$Nsae <- NULL


##     setkeyv(shc, c('nbatch', 'simNum', whichDo))

##     shc <- merge(parms, shc, all.y=T, by = c('nbatch'))
##     setkeyv(shc, c('nbatch', 'simNum', whichDo))

##     ## shc[,.N, list(propInTrial, ord, gs, trialStartDate, trial, simNum, nbatch, avHaz, vaccDay_noEV)] # number of hazard categories
##     ## shc[,.N, list(propInTrial, ord, gs, trialStartDate, trial, simNum, nbatch, avHaz, vaccDay_noEV, ihaz0Cat)] # down to unique

##     ## shc[, .N, list(propInTrial, trialStartDate, simNum, avHaz, vaccDay_noEV, ihaz0Cat)] ## 7 simulation types
##     ## unique(parms[propInTrial==.025 & trialStartDate=='2014-10-01' & avHaz==''][,2:12,with=F])

##     browser()

    
##     ## infection risk within strata
## shc[, risk_EV:=Ninf_EV/N]
## shc[, risk_noEV:=Ninf_noEV/N]

##     ## for averted: NT risk comparator is same risk group & given they're never vaccinated
## shc[, risk_NT:= risk_noEV[trial=='NT' & vaccDay_noEV==Inf], list(propInTrial, trialStartDate, simNum, avHaz, ihaz0Cat)]
##     ## for spent: VR risk comparator should be within same risk group marginal over vaccination group, vs that risk group marginal over vaccination group again. Cannot look at comparisons 
## shc[, risk_VR:= risk_noEV[trial=='VR'], list(propInTrial, trialStartDate, simNum, avHaz, vaccDay_noEV, ihaz0Cat)]
##     ## figuring out e* metrics here

##         shc[N>0]
    
##     shc[, N_NT:= Ninf_noEV[trial=='NT'], list(propInTrial, trialStartDate, simNum, avHaz, ihaz0Cat)]


##     shc[, N_VR:= Ninf_noEV[trial=='VR'], list(propInTrial, trialStartDate, simNum, avHaz, vaccDay_noEV, ihaz0Cat)]
##     ## Averted infections *marginal on randomization assignment*
##     shc[, infAvert_EV := sum(N_NT) - sum(Ninf_EV), list(propInTrial, trialStartDate, simNum, avHaz, ihaz0Cat)]
##     shc[, infAvert_noEV := sum(N_NT) - sum(Ninf_noEV), list(propInTrial, trialStartDate, simNum, avHaz, ihaz0Cat)]
##     shc[, infAvertPC_EV := infAvert_EV/sum(N), list(propInTrial, trialStartDate, simNum, avHaz, ihaz0Cat)]
##     shc[, infAvertPC_noEV := infAvert_noEV/N]
##     shc[is.nan(infAvertPC_EV), infAvertPC_EV := 0]
##     shc[is.nan(infAvertPC_noEV),infAvertPC_noEV := 0]
##     ## shc[, infAvertProp_EV := infAvert_EV /N_NT]
##     ## shc[, infAvertProp_noEV := infAvert_noEV /N_NT]
##     ## Spent infections *conditional* on randomization assignment
##     shc[, infSpent_EV := Ninf_EV - N_VR]
##     shc[, infSpent_noEV := Ninf_noEV - N_VR]
##     shc[, infSpentPC_EV := infSpent_EV/N]
##     shc[, infSpentPC_noEV := infSpent_noEV/N]
##     shc[is.nan(infSpentPC_EV), infSpentPC_EV := 0]
##     shc[is.nan(infSpentPC_noEV),infSpentPC_noEV := 0]
##     ## shc[, infSpentProp_EV := infSpent_EV /N_NT]
##     ## shc[, infSpentProp_noEV := infSpent_noEV /N_NT]

##     shc[, lab:=trial]
##     shc[trial=='RCT' & gs==T, lab:=paste0(lab,'-gs')]
##     shc[trial=='RCT' & ord=='TU', lab:=paste0(lab,'-rp')]
##     shc[,lab:=as.factor(lab)]

##     ## **Eventually want fraction of infections that were used in analysis? not sure think about
##     ## **Careful because now splitting the information between randomized vaccine dates, need to add
##     ## it back up for something meaningful (it's vacc & cont information together that provides any
##     ## information really)

##     shc[, NinfFrac_EV := Ninf_EV/sum(Ninf_EV), list(propInTrial, trialStartDate, simNum, nbatch, lab, avHaz)]
##     shc[, NinfFrac_noEV := Ninf_noEV/sum(Ninf_noEV), list(propInTrial, trialStartDate, simNum, nbatch, lab, avHaz)]
    
##     infAvertPow <- shc[trial!='NT', list(.N
##                          , infAvert_EV = mean(infAvert_EV), infAvert_noEV = mean(infAvert_noEV)
##                          , infAvertPC_EV = mean(infAvertPC_EV), infAvertPC_noEV = mean(infAvertPC_noEV)
##                          , infSpent_EV = mean(infSpent_EV), infSpent_noEV = mean(infSpent_noEV)
##                          , infSpentPC_EV = mean(infSpentPC_EV), infSpentPC_noEV = mean(infSpentPC_noEV) 
##                          , NinfFrac_EV = mean(NinfFrac_EV), NinfFrac_noEV = mean(NinfFrac_noEV)),
##                        ##infAvertProp = mean(infAvertProp),  infAvert_noEVProp = mean(infAvert_noEVProp)
##                        list(propInTrial, trialStartDate, lab, avHaz, vaccDay_noEV, ihaz0Cat)]

##     infAvertPow <- merge(infAvertPow
##                        , finit[cat=='allFinalEV' & trial!='NT', list(powerEff=mean(vaccGood[vaccEff>0]), tcalEff=mean(tcal[vaccEff>0])), 
##                                list(propInTrial, trialStartDate, lab, avHaz)]
##                        , by = c('propInTrial', 'trialStartDate', 'lab', 'avHaz'))
##     infAvertPow[lab=='VR', powerEff:=0]

##     infAvertPow[, powEffFrac_EV:= NinfFrac_EV * powerEff]
##     infAvertPow[, powEffFrac_noEV:= NinfFrac_noEV * powerEff]
##     infAvertPow[lab=='VR', c('powEffFrac_EV', 'powEffFrac_noEV') := 0]
##     infAvertPow[, powEff_per_infSpent_EV:=powEffFrac_EV/infSpent_EV]
##     infAvertPow[, powEff_per_infSpent_noEV:=powEffFrac_noEV/infSpent_noEV]
## })


makePowTable <- function(finTrials, verbose=0, doSave=T) {
    if(verbose==1) browser()
    pf <- finTrials[, list(
        nsim = length(stopped)
      , stopped = mean(stopped)
      , vaccGood = mean(vaccGood)
      , stopDay = mean(tcal)
      , caseTot = mean(vaccCases+contCases)
      , caseC = mean(contCases)
      , caseV = mean(vaccCases)
      , cvr = mean(cvr)
      , cvrNAR = mean(cvr, na.rm=T)
      , mean = mean(mean)
      , meanNAR = mean(mean, na.rm=T)
      , vaccBad = mean(vaccBad)
      , stoppedNAR = mean(stopped,na.rm=T)
      , vaccGoodNAR = mean(vaccGood,na.rm=T)
      , vaccBadNAR = mean(vaccBad,na.rm=T)
      , PHUNAR = mean(PHU, na.rm=T)
      , meanErr = mean(err)
      , meanBump = mean(bump)
        ## , caseTot = mean(caseCXimmGrpEnd + caseVXimmGrpEnd )
        ## , caseC = mean(caseCXimmGrpEnd)
        ## , caseV = mean(caseVXimmGrpEnd)
    ),
                    list(vaccEff, trial, gs, propInTrial, ord, delayUnit, mod, immunoDelay,trialStartDate
                       , weeklyDecay, cvWeeklyDecay, cvClus) #, cvClusTime)
                    ]


    if(verbose==2) browser()

    coerce.numeric <- function(dat = pf, variables = c('propInTrial', 'delayUnit','immunoDelay')) within(dat, {
        for(vv in variables) {
            if(dat[,class(get(vv))]=='character') dat[,vv:= as.numeric(get(vv)), with=F]
            if(dat[,class(get(vv))]=='factor') dat[,vv:= as.numeric(levels(get(vv))[get(vv)]), with = F]
        }
        rm(vv)
    })
    coerce.numeric(pf)

    pf[,trial:=factor(trial)]
    ## Get bias from means done on a PHU scale
    pf$RHxPHUNAR <- pf[, -PHUNAR/(PHUNAR-1)]
    pf$meanXPHUNAR <- pf[, 1 - RHxPHUNAR]
    pf$biasNAR <- pf[, meanXPHUNAR - as.numeric(vaccEff)]

    ## Formatting stuff
    front <- c('mod','vaccEff','stoppedNAR','vaccGoodNAR','cvrNAR','biasNAR',
               'nsim','meanErr','propInTrial','vaccBad','cvr','stopped','vaccGood')
    setcolorder(pf, c(front, setdiff(names(pf), front)))
    pf <- pf[!(trial=='FRCT' & delayUnit==0) & !(ord=='TU' & delayUnit==0)] ## redundant
    pf$trialStartDate <- as.Date(pf$trialStartDate)
    pf[mod=='coxME', mod:='CoxME']
    pf$mod <- factor(pf$mod, levels=unique(pf$mod))
    pf$order <- pf$ord
    pf[delayUnit==0, order:='simultaneous instant']
    pf$design <- pf$trial
    levels(pf$design)[levels(pf$design) == 'SWCT'] <- 'SWT'
    levels(pf$order)[2] <- 'time-updated'

    pf[, pit:=factor(paste0(propInTrial*100,'%'))]
    pf[, pit:=factor(pit, levels = c('2.5%','5%','7.5%','10%'), ordered = T)]
    baseMods <- c('Cox PH Frailty'
                , 'Poisson GLM\n no cluster effects'
                , 'Poisson GLM \nwith fixed effects by cluster')
    pf$model <- pf$mod
    levels(pf$model) <- paste0(rep(c('', 'bootstrap over\n', 'permutation test over\n'),each=3), rep(baseMods,3))
    if(doSave)
        save(pf, file=file.path('Results',paste0('pf_',thing,'.Rdata')))
    return(pf)
}
