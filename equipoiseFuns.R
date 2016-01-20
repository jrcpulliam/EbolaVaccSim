## Uncertainty about vaccine effectiveness

equipFxn <- function(N=1
                     , seed=1
                     , probVaccWorks = .7
                     , EfficacyRNG = function(N) runif(N, .5, .9) ## pdf of efficacy give any efficacy
                     , saeRNG = function(N) runif(N, 10^-5, 5*10^-4) ## serious adverse effect
                     , probVaccIncreasesCFR = 0 ## probablity a vaccine is not efficacious AND increases the case fatality rate
                     , badVaccCFR_OR = 2 ## odds ratio of dying given EVD if vaccinted versus if not, if above true                     
                     ) {
    vaccWorks <- rbinom(N, 1, probVaccWorks)
    vaccEff <-  EfficacyRNG(N) * vaccWorks
    vaccBad <- rbinom(N, 1, probVaccIncreasesCFR / (1-probVaccWorks)) * (1 - vaccWorks) ## if vaccine doesn't work
    pSAE <- saeRNG(N)
    cfrOR <- badVaccCFR_OR * vaccBad
    return(data.table(vaccEff, pSAE, cfrOR))
}

## set.seed(1)
## vaccProp1 <- equipFxn(5000)
## save(vaccProp1, file='data/vaccProp1.Rdata')

## pdf('Figures/EquipPrior.pdf')
## par(mfrow=c(2,1))
## hist(vaccProp1$vaccEff, xlab = 'vaccine efficacy', col = 'black', bty = 'n')
## hist(vaccProp1$pSAE, xlab = 'pSAE', col = 'black', bty = 'n')
## dev.off()

## equipFxn(10)

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
