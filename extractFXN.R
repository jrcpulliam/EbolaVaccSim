popNms <- c("Spop", "SpopH", "SpopEvents", "SEVpopEvents")
dparms0 <- c('trial','gs','doSL','propInTrial','nbsize', 'avHaz', 'clusSize', 'numClus', ## **automate 150,300,6000** change later on new batches
             'ord','reordLag','delayUnit' ,'immunoDelay','trialStartDate', 'HazTrajSeed'                          
             )
moveFront <- function(dt, fnames) setcolorder(dt, c(fnames, names(dt)[!names(dt) %in% fnames]))

quantcut <- function(x, qs = seq(0,1, l = 6)) {
    brks <- unique(quantile(x, qs))
    if(length(brks)==1) brks <- brks*c(.99,1/.99) ## in case it's one number, just make an interval that includes it
    return(as.numeric(cut(x, brks, include.lowest = T)))
}

procAll <- function(tidDo, verbose=0) {
    if(verbose>0) browser()
    nbtd <- parmsMat[tid==tidDo,rcmdbatch]
    resList <- extractSims(thing, verb=0, maxbatches=NA, nbatchDo=nbtd, indivLev = T, mc.cores=48)
    resList <- procResList(resList, verb=0)
    resList <- procMetaParms(resList)
    resList <- procIrskSpent(resList)
    save(resList, file=file.path('BigResults',paste0(thing, '-', tid, '.Rdata')))
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
                setkey(SpopH, simNum, Oc)
                Spop <- data.table(nbatch = nbatch, Spop)
                SpopH <- data.table(nbatch = nbatch, SpopH)
            })
            essence <- c(essence, riskStratList)
        }
        return(essence)
    }
}

## Extract from simulations from multiple cores
extractSims <- function(thing
                      , dparms = dparms0
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
        fls <- fls[fns %in% nbatchDo]
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
    return(resList)
}

procResList <- function(resList, verbose = 0) {
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
    return(resList)
}

## ####################################################################################################
## ## spent vs power
## ptab <- merge(finTrials[,!"sim", with=F], parms, by = 'nbatch')[, list(pid, lab, simNum, vaccEff, mean, lci, uci, vaccGood, vaccBad, stopped, cvr)]
## triSumm <- merge(SpopWorst, ptab[vaccEff>0 & pid < 6, list(power = mean(vaccGood)), list(pid)], by = 'pid')
## triSumm <- merge(triSumm, Spop[pid<6,list(totCase = sum(infectDay_EV<336), .N), list(pid, simNum)][,list(totCase = mean(totCase)),pid], by = 'pid')

## add risk spent & averted for individuals to Spop & also create punq & parms
procMetaParms <- function(resList)  within(resList, {
    ## unique parameters
    punq <- unique(parms[,-1,with=F])
    punq <- data.table(pid = 1:nrow(punq), punq)
    punq[, lab:=trial]
    punq[trial=='RCT' & gs==T, lab:=paste0(lab,'-gs')]
    punq[trial=='RCT' & ord=='TU', lab:=paste0(lab,'-rp')]
    punq[,lab:=as.factor(lab)]
    punq[,lab:=factor(lab, levels =levels(lab)[c(1:3,5,4,6:7)])]
    setkey(punq, pid)
    ## create parms
    parms <- merge(punq, parms, by = names(punq)[!names(punq) %in% c('lab','pid')])
    setkey(parms, pid, nbatch)
    setcolorder(parms, c('pid','nbatch', names(parms)[!names(parms) %in% c('pid','nbatch')]))
    Spop <- merge(parms[,list(pid, lab, nbatch)], Spop, by = 'nbatch')
    setcolorder(Spop, c('pid','lab','nbatch', names(Spop)[!names(Spop) %in% c('pid','lab', 'nbatch')]))
    Spop[,arm:=c('vacc','cont')[as.numeric(vaccDay==Inf)+1]]
    setkey(Spop, Oi, pid, nbatch, simNum)
    ## Spop[,quantile(indivRR, c(.025,.975))]
    Spop[,cumRisk:=1-exp(-cumHaz)] 
    Spop[,cumRisk_EV:=1-exp(-cumHaz_EV)]
    ## Calculate risk spent so that can tabulate expected # of people spending a certain amount of risk in a given simulation
    Spop[, spent:= cumRisk - cumRisk[lab=='VR'],list(simNum, Oi)]
    Spop[, avert:= cumRisk[lab=='NT'] - cumRisk,list(simNum, Oi)]
    Spop[, spent_EV:= cumRisk_EV - cumRisk[lab=='VR'],list(simNum, Oi)]
    Spop[, avert_EV:= cumRisk[lab=='NT'] - cumRisk_EV,list(simNum, Oi)]
})

procExpRiskSpent <- function(resList, threshold = .2, breaks = seq(-.5, 1, by = .01)) within(resList, {
    ## frequency distribution of risk spent for each simulation, then will get expected # of people within each risk spent bin across 2000 scenarios
    HpopInd <- Spop[!lab %in% c('VR','NT'), hist(spent, breaks=breaks, plot=F)[c('mids','counts')], list(pid, simNum)]
    HpopInd_EV <- Spop[!lab %in% c('VR','NT'), hist(spent_EV, breaks=breaks, plot=F)[c('mids','counts')], list(pid, simNum)]
    HpopInd <- merge(HpopInd, HpopInd_EV, by = c('pid','simNum','mids'), suffixes=c('','_EV'))
    rm(HpopInd_EV)
    Hpop <- HpopInd[, list(counts=median(counts), counts_EV=median(counts_EV)), list(pid, mids)]
    Hpop <- merge(Hpop, punq, by = 'pid')
    ## Find the simulation with the most # of people above risk spent threshold
    SpopThresh <- Spop[!lab %in% c('VR','NT'), list(above = sum(spent > threshold), above_EV = sum(spent_EV > threshold), .N), list(pid, simNum)]
    SpopWorst <- merge(SpopThresh[, list(above = max(above), above_EV = max(above_EV), .N), list(pid)], punq, by = 'pid')
})

## **automate 150,300,6000** change later on new batches
procIrskSpent <- function(resList) within(resList, {
    ## table by individual the average infection risk in each design
    irskMarg <- Spop[, list(.N, inf = mean(cumRisk), inf_EV = mean(cumRisk_EV) 
                         ,  indivRR=unique(indivRR), Oc=Oc[1], type = 'marg', arm = NA), 
                     list(pid,Oi)]
    irskCond <- Spop[!lab %in% c('VR','NT','SWCT'), ## conditional on control/vacc randomization assignment
                     list(.N, inf = mean(cumRisk), inf_EV = mean(cumRisk_EV)
                        , indivRR=unique(indivRR), Oc=Oc[1], type = 'cond'), 
                     list(pid,Oi,arm)]
    ## maximum spent given worst possible vaccination order (only for random ordered trials)
    irskMax <- Spop[!lab %in% c('VR','NT') & arm=='vacc' & !grepl('rp',lab) & vaccDay==max(vaccDay[arm=='vacc']),
                    list(.N, inf = mean(cumRisk), inf_EV = mean(cumRisk_EV), indivRR=unique(indivRR), Oc=Oc[1], arm='vacc', type='max'),
                    list(pid,Oi,vaccDay)]
    ## conditional on exact vaccination day (note borrowing control info across all vaccDay==Inf, i.e. whether the cluster they're in is vaccinated early/late)
    irskCondvd <- Spop[!lab %in% c('VR','NT'),
                       list(.N, inf = mean(cumRisk), inf_EV = mean(cumRisk_EV), indivRR=unique(indivRR), Oc=Oc[1], type='condvd'),
                       list(pid,Oi,arm,vaccDay)]
    irsk <- rbindlist(list(irskMarg, irskCond, irskMax, irskCondvd), use.names=T, fill=T)
    rm(irskMarg, irskCond, irskMax, irskCondvd)
    setkey(irsk, Oi, pid)
    irsk <- merge(punq[,list(pid,gs, lab)], irsk, by = 'pid')
    moveFront(irsk, c('pid','lab','gs','type','arm','vaccDay'))
    ## irsk
    ## irsk[, unique(type), lab]
    ## irsk[Oi==1]
    ## Average risk spent versus marginal VR 
    irsk[, spent   := inf    - inf[lab=='VR'], list(Oi)]
    irsk[, spent_EV:= inf_EV - inf[lab=='VR'], list(Oi)]
    ## Averted risk versus marginal NT
    irsk[, avert   :=    inf[lab=='NT'] - inf, list(Oi)]
    irsk[, avert_EV:= inf[lab=='NT'] - inf_EV, list(Oi)]
    setkey(irsk, Oi, pid)
    ## ########
    ## Examine results
    ## irsk[type!='condvd', list(inf = mean(inf), inf_EV = mean(inf_EV)), list(lab, arm, type)][order(lab)] ## mean infection rate without vaccine
    ## sptmp <- irsk[type!='condvd', list(spent = mean(spent), spent_EV = mean(spent_EV)), list(lab, arm, type)][order(lab)] #
    ## for(ii in 4:5) sptmp[[ii]] <- formatC(100*signif(sptmp[[ii]],2)) ## % risk spent on average
    ## sptmp ## % risk (multiplied by 100)
    ## ##############################################################################
    ## get clusters in risk-order
    cOrd <- irsk[lab=='NT', list(inf=mean(inf)), Oc][order(-inf)]
    cOrd[,cluster:=1:nrow(cOrd)]
    irsk$cluster <- NULL
    irsk <- merge(irsk, cOrd[,list(Oc, cluster)], by = 'Oc')
    irsk[,cluster:=factor(cluster)]
    ## order individuals within clusters by risk for ease of display
    iord <- irsk[lab=='NT',list(Oi,cluster,inf)][order(cluster,-inf)]
    iord[,ordShow:=1:6000]
    irsk$ordShow <- NULL
    irsk <- merge(irsk, iord[,list(Oi,ordShow)], by = 'Oi')
    ## what arm to label individuals as for RCTs?
    irsk[, armShown:=arm]
    ## if(length(unique(parms[,clusSize]))>1) stop("more than 1 clusSize value in simulation batch")
    ## clusSize <- parms[1,clusSize]
    clusSize <- 300
    ## if(length(unique(parms[,numClus]))>1) stop("more than 1 numClus value in simulation batch")
    ## numClus <- parms[1,numClus]
    numClus <- 20
    irsk[grepl('RCT',lab) & as.numeric((Oi-1) %% (clusSize) < clusSize/2), armShown:='cont']
    irsk[grepl('RCT',lab) & as.numeric((Oi-1) %% (clusSize) >= clusSize/2), armShown:='vacc']
    ## order individuals within clusters & armShown by risk for ease of display
    iord <- irsk[lab=='NT',list(Oi,cluster,inf)]
    iord <- merge(iord, unique(irsk[lab=='RCT' & arm==armShown, list(Oi, armShown)]))
    iord <- iord[order(cluster,armShown, -inf)]
    iord[,ordShowArm:=1:(clusSize*numClus)]
    irsk$ordShowArm <- NULL
    irsk <- merge(irsk, iord[,list(Oi,ordShowArm)], by = 'Oi')
    irsk[lab=='SWCT', ordShowArm:=ordShow]
    irsk <- irsk[order(Oc,indivRR,lab)]
    irsk[type=='cond' &  !lab %in% c('VR','NT') & pid==2 & Oi==1]
    ## need exemplar SWCT: pick arbitrary example of cluster-day assignment to display
    ## **(later do probably do exemplars for all of them with exact same randomization throughout)
    clusVD <- irsk[(arm==armShown & type=='condvd' &  grepl('SWCT',lab)),list(Oc = sample(1:numClus,numClus), vaccDay=unique(vaccDay))]
    setkey(clusVD, Oc, vaccDay)
    irsk[,exmpl:=F]
    setkey(irsk, Oc, vaccDay)
    irsk[clusVD][type=='condvd' & 'SWCT'==lab][["exmpl"]] <- rep(T,6000)
    irsk[clusVD][type=='condvd' & 'SWCT'==lab]
    setkey(irsk, Oi, pid)
    rm(iord, cOrd)
})

## Make a data table of risk-strata level info on infections spent/averted & power contributed by
## that group within the trial
makeInfPow <- function(resList, doSave=T, verbose = 0) within(resList, {
    if(verbose>0) browser()

    ## unique identifier for trial population & for parameter set
                   
    
    ##     if(doSave) {
    ##         save(resList, file=file.path('BigResults', paste0(thing, '.Rdata')))
    ##     }

})


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
