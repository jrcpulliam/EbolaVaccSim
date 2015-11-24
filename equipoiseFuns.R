## Uncertainty about vaccine effectiveness

equipFxn <- function(N=1
                     , seed=1
                     , probVaccWorks = .7
                     , efficacyRange = c(.5,.9)
                     , SAERange = c(10^-5, 10^-4) ## serious advers effect
                     , probVaccIncreasesCFR = .1 ## probablity a vaccine is not efficacious AND increases the case fatality rate
                     , badVaccCFR_OR = 2 ## odds ratio of dying given EVD if vaccinted versus if not, if above true                     
                     ) {
    vaccWorks <- rbinom(N, 1, probVaccWorks)
    vaccEff <- runif(N, efficacyRange[1], efficacyRange[2]) * vaccWorks
    vaccBad <- rbinom(N, 1, probVaccIncreasesCFR / (1-probVaccWorks)) * (1 - vaccWorks)
    SAE <- runif(N, SAERange[1], SAERange[2])
    cfrOR <- badVaccCFR_OR * vaccBad
    return(data.table(vaccEff, SAE, cfrOR))
}

## equipFxn(100)

equiCalc <- function(parms) within(parms, {
    browser()

    ls()
    ## Difference between trial & TU vacc rollout
    finInfo[cat=='allFinal_VR', caseTot] - finInfo[cat=='allFinal_EV', caseTot]

})

####################################################################################################
## Build counterfactual simulation sets for comparison with factual simulations
## Important to make sure that they match the seeds (can use indivHaz at end to do a check)
cfSims <- function(parms) {
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
    return(parms)
}
## parms$VRpopH[idByClus %in% c(1,151)]

simN_CFs <- function(seed = 1, parms=makeParms(), N = 2, verbFreq=10) {
    set.seed(seed)
    finInfo <- data.frame(NULL)
    for(ss in 1:N) {
        if(parms$verbose>0 & (ss %% verbFreq == 0)) print(paste('on',ss,'of',N))
        if(parms$verbose>.5 & (ss %% 1 == 0)) print(paste('on',ss,'of',N))
        if(parms$verbose==2) browser()
        res <- simTrial(parms)
        res <- cfSims(res)
        res$doCFs <- T
        res <- finInfoFxn(res)
        ## compile results from the final interim analysis (or all statistical analyses for a single fixed design)
        finITmp <- data.table(sim = ss, res$finInfo)
        finInfo <- rbind(finInfo, finITmp)
        ## res <- equiCalc(res)
        rm(res)
        gc()
    }
        return(list(finInfo=finInfo))
}

## simN_CFs(1, makeParms(verbose=1))
