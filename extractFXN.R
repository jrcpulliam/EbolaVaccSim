popNms <- c("Spop", "SpopH", "SpopEvents", "SEVpopEvents")
dparms0 <- c('trial','gs','doSL','propInTrial','nbsize', 'avHaz',
                            'ord','reordLag','delayUnit' ,'immunoDelay','trialStartDate', 'HazTrajSeed'                          
                                   )

quantcut <- function(x, qs = seq(0,1, l = 6)) {
    brks <- unique(quantile(x, qs))
    if(length(brks)==1) brks <- brks*c(.99,1/.99) ## in case it's one number, just make an interval that includes it
    return(as.numeric(cut(x, brks, include.lowest = T)))
}

extractOneSim <- function(fileNm
                        , dparms = dparms0
                        , verbose = 0
                        , verbFreq = 100
                        , indivLev = F, hazBrks = 10^seq(-12,1, by = 1)
                          ) {
    load(fileNm)
    riskStratList <-  finInfoList <- finModList <- parmsList <- list(NULL)
    if(exists('sim')) {
        if(verbose==1) browser()
        nbatch <- as.numeric(gsub("[^0-9]", "", fileNm))
        if(nbatch %% verbFreq == 0) print(nbatch)
        sim$parms[['trialStartDate']] <- as.character(sim$parms[['trialStartDate']])
        parmsList <- data.table(nbatch = nbatch, t(unlist(sim$parms[dparms])))
        finModList <- data.table(nbatch = nbatch, sim$sim$finMods)
        finInfoList <- data.table(nbatch = nbatch, sim$sim$finInfo)
        if(indivLev) {
            riskStratList <- with(sim$sim$Spops, {
                setkey(Spop, simNum, indiv, cluster)
                setkey(SpopH, simNum, cluster)
                setkey(SpopEvents, simNum, indiv)
                notCF <- nrow(SEVpopEvents)>0
                if(notCF) setkey(SEVpopEvents, simNum, indiv)
                ## get first hazard & indivRR for every infected individual
                unqvars <- c('simNum')
                SpopH[, avHaz:= mean(clusHaz), list(simNum, cluster)]
                SpopH[, haz0:= clusHaz[day==0], list(simNum, cluster)]
                Spop <- merge(Spop, SpopH[day==0, list(simNum, cluster, haz0, avHaz)], by = c(unqvars, 'cluster'))
                Spop[,haz0Q:=quantcut(haz0), unqvars] ## quantiles *within* simulations
                Spop[,ihaz0Q:=quantcut(haz0*indivRR), unqvars]
                Spop[,haz0Cat:=cut(haz0, hazBrks, include.lowest = T)] ## absolute breaks in hazard
                Spop[,ihaz0Cat:=cut(haz0*indivRR, hazBrks, include.lowest = T)]
                ## number infected & SAE overall, with individual level data tracked
                trackUntilDay <- sim$sim$finInfo[cat=='allFinalEV',atDay][1]

                if(notCF) {

                    infPop_EV <- merge(SEVpopEvents[, list(simNum, indiv, infectDay)], Spop, by=c(unqvars,'indiv'))
                    saePop_EV <- merge(SEVpopEvents[, list(simNum, indiv, SAE)], Spop, by=c(unqvars,'indiv'))
                }

                infPop_noEV <- merge(SpopEvents[, list(simNum, indiv, infectDay)], Spop, by=c(unqvars,'indiv'))
                saePop_noEV <- merge(SpopEvents[, list(simNum, indiv, SAE)], Spop, by=c(unqvars,'indiv'))


                ## ##################################################
                ## To check that individual-level data match pop-aggregate compare infPop & finInfo
                ## vtmp <- merge(infPop[,list(caseTot = .N), unqvars], 
                ##               sim$sim$finInfo[cat=='allFinalEV', list(simNum, caseTot)], by=unqvars)
                ## vtmp[, range(caseTot.x-caseTot.y)] ## should match
                ## vtmp <- merge(infPop_noEV[,list(caseTot = .N), unqvars], 
                ##               sim$sim$finInfo[cat=='allFinal_noEV', list(simNum, caseTot)], by=unqvars)
                ## vtmp[, range(caseTot.x-caseTot.y)] ## should match
                ## ##################################################
                ## table the number of infections by each type
                vars <- c('haz0Q', 'haz0Cat', 'ihaz0Cat')
                riskStratList <- list()
                for(vv in vars) {
                    if(notCF) {
                        ## infection tallies by risk level at beginning of trial
                        itmp <- merge(infPop_EV[,list(.N, Ninf=length(simNum[infectDay<Inf])), c(unqvars,vv)] 
                                    , infPop_noEV[,list(Ninf=length(simNum[infectDay<Inf])), c(unqvars,vv)] 
                                    , by = c('simNum',vv), all=T, suffixes = c('_EV','_noEV'))
                        stmp <- merge(saePop_EV[,list(Nsae=length(simNum[SAE==1])), c(unqvars,vv)] 
                                    , saePop_noEV[,list(Nsae=length(simNum[SAE==1])), c(unqvars,vv)] 
                                    , by = c('simNum',vv), all=T, suffixes = c('_EV','_noEV'))
                    }else{
                        itmp <- infPop_noEV[,list(.N, Ninf=length(simNum[infectDay<Inf])), c(unqvars,vv)] 
                        stmp <- saePop_noEV[,list(Nsae=length(simNum[SAE==1])), c(unqvars,vv)] 
                    }
                    tmp <- merge(itmp, stmp, by = c('simNum',vv), suffixes = c('inf','sae'), all=T)
                    Nnms <- colnames(tmp)[grepl('N_',colnames(tmp))]
                    for (col in Nnms) set(tmp, which(is.na(tmp[[col]])), col, 0)
                    tmp$nbatch <- nbatch
                    setcolorder(tmp, c('nbatch', names(tmp)[names(tmp)!='nbatch']))
                    riskStratList[[paste0('strat_',vv)]] <- tmp
                }
                return(riskStratList)
            } ) ## with(Spops)
        } ## indivLevrisk
    }
    essence <- list(parms=parmsList, finTrials=finModList, finInfo=finInfoList)
    essence <- c(essence, riskStratList)
    ## sim only used for managing within core-batches so delete it
    for(vv in names(essence)) if('sim' %in% names(essence[[vv]])) essence[[vv]]$sim <- NULL 
    return(essence)
}

## Extract from simulations from multiple cores
extractSims <- function(thing
                      , dparms = dparms0
                      , doSave=T
                      , verbose = 0
                      , maxbatches = NA
                      , mc.cores = 48
                      , indivLev = F, hazBrks = 10^c(-12:0)
                        ) {
    if(verbose==1) browser()
    batchdirnm <- file.path('BigResults',thing)
    fls <- list.files(batchdirnm, pattern=thing, full.names = T)
    nbatch <- length(fls)
    if(!is.na(maxbatches)) {
        nbatch <- maxbatches
        fls <- fls[1:nbatch]
    }
    print(paste0('extracting from ', nbatch, ' files'))
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

        finit[order(gs), list(tcalMean = mean(tcal), power = mean(vaccGood), length(tcal)),
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
