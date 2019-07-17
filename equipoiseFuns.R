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
cfSims <- function(parms, seed=1) {
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
    parms$InfTimesLs <- list()
    for(cf in c('NTpop','VRpop')) {
        ## CONTROL SEEDS HERE
        parms$InfTimesLs[[cf]] <- list()
        length(parms$InfTimesLs[[cf]]) <- parms$numCFs
        for(cc in 1:parms$numCFs) {
            parms$InfTimesLs[[cf]][[cc]] <- data.table(seed = seed, cf=cf, cc=cc, simInfection(parms, cfNum=cc)$indivInfDays) ## simulate infection
        }
        parms$InfTimesLs[[cf]] <- rbindlist(parms$InfTimesLs[[cf]])
    }
    parms$InfTimesLs <- rbindlist(parms$InfTimesLs)
    return(parms)
}
## parms$VRpopH[idByClus %in% c(1,151)]

## Compress infection time information since we may want to track how
compInfTimes <- function(parms, whichDo = c('NTpop','VRpop')) within(parms, {
    InfTimes <- list()
    for(ww in whichDo) {
        tmp <- get(ww)
        InfTimes[[ww]] <- tmp[infectDay!=Inf, list(indiv, infectDay, indivRR)]
    }
    rm(tmp,ww)
})

## Need distinct set of counterfactuals for these parameters only (which determine population,
## hazard, & maximum rollout speed ). Others excluded here only affect trial design.
CFparms <- c('trialStartDate'
           , 'numClus'
           , 'clusSize'
           , 'delayUnit'
           , 'hazType'
           , 'nbsize'
           , 'propInTrial'
           , 'mu'
           , 'cvClus'
           , 'cvClusTime'
           , 'sdLogIndiv'
           , 'vaccEff'
           , 'trackUntilDay'
           , 'immunoDelay'
           , 'weeklyDecay'
           , 'cvWeeklyDecay'
           , 'hazIntUnit'
           , 'reordLag')
