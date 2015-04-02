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

equipFxn(100)
