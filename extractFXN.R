popNms <- c("Spop", "SpopH", "SpopEvents", "SEVpopEvents")
dparms0 <- c('trial','gs','doSL','propInTrial','nbsize', 
                            'ord','reordLag','delayUnit' ,'immunoDelay','trialStartDate'
                          , 'weeklyDecay', 'cvWeeklyDecay', 'cvClus', 'cvClusTime', 'avHaz'
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
                        , indivLev = F, hazBrks = 10^c(-12:0)
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
                setkey(Spop, sim, simNum, indiv, cluster)
                setkey(SpopH, sim, simNum, cluster)
                setkey(SpopEvents, sim, simNum, indiv)
                notCF <- nrow(SEVpopEvents)>0
                if(notCF) setkey(SEVpopEvents, sim, simNum, indiv)
                ## get first hazard & indivRR for every infected individual
                unqvars <- c('sim', 'simNum')
                SpopH[, avHaz:= mean(clusHaz), list(sim, simNum, cluster)]
                SpopH[, haz0:= clusHaz[day==0], list(sim, simNum, cluster)]
                Spop <- merge(Spop, SpopH[day==0, list(sim, simNum, cluster, haz0, avHaz)], by = c(unqvars, 'cluster'))
                Spop[,haz0Q:=quantcut(haz0), unqvars] ## quantiles *within* simulations
                Spop[,ihaz0Q:=quantcut(haz0*indivRR), unqvars]
                Spop[,haz0Cat:=cut(haz0, hazBrks, include.lowest = T)] ## absolute breaks in hazard
                Spop[,ihaz0Cat:=cut(haz0*indivRR, hazBrks, include.lowest = T)]
                ## number infected & SAE overall, with individual level data tracked
                trackUntilDay <- sim$sim$finInfo[cat=='allFinalEV',atDay][1]
                if(notCF) {
                    infPop <- merge(SEVpopEvents[infectDay<trackUntilDay, list(sim, simNum, indiv)], Spop, by=c(unqvars,'indiv'))
                    saePop <- merge(SEVpopEvents[SAE>0, list(sim, simNum, indiv)], Spop, by=c(unqvars,'indiv'))
                }
                infPop_noEV <- merge(SpopEvents[infectDay<trackUntilDay, list(sim, simNum, indiv)], Spop, by=c(unqvars,'indiv'))
                saePop_noEV <- merge(SpopEvents[SAE>0, list(sim, simNum, indiv)], Spop, by=c(unqvars,'indiv'))
                ## ##################################################
                ## To check that individual-level data match pop-aggregate compare infPop & finInfo
                ## vtmp <- merge(infPop[,list(caseTot = .N), unqvars], 
                ##               sim$sim$finInfo[cat=='allFinalEV', list(sim, simNum, caseTot)], by=unqvars)
                ## vtmp[, range(caseTot.x-caseTot.y)] ## should match
                ## vtmp <- merge(infPop_noEV[,list(caseTot = .N), unqvars], 
                ##               sim$sim$finInfo[cat=='allFinal_noEV', list(sim, simNum, caseTot)], by=unqvars)
                ## vtmp[, range(caseTot.x-caseTot.y)] ## should match
                ## ##################################################
                ## table the number of infections by each type
                vars <- c('haz0Q', 'haz0Cat', 'ihaz0Cat')
                riskStratList <- list()
                for(vv in vars) {
                    if(notCF) {
                        itmp <- merge(infPop_noEV[,.N, c(unqvars,vv)] ## infection tallies by risk level
                                    , infPop[,.N, c(unqvars,vv)]
                                    , by = c('sim','simNum',vv), all=T, suffixes = c('_EV','_noEV'))
                        stmp <- merge( saePop_noEV[,.N, c(unqvars,vv)] ## sae tallies by risk level
                                    , saePop[,.N, c(unqvars,vv)]
                                    , by = c('sim','simNum',vv), all=T, suffixes = c('_EV','_noEV'))
                    }else{
                        itmp <- infPop_noEV[,list(N_noEV = .N), c(unqvars,vv)]
                        stmp <- saePop_noEV[,list(N_noEV = .N), c(unqvars,vv)]
                    }
                    tmp <- merge(itmp, stmp, by = c('sim','simNum',vv), suffixes = c('inf','sae'))
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
    fls <- list.files(batchdirnm, pattern='.Rdata', full.names = T)
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
        if(vv!='parms') setkey(resList[[vv]], nbatch ,sim, simNum)
    }

    resList <- within(resList, {
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
