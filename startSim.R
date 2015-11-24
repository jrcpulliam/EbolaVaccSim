if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R', 'MungeFuns.R','CoxFxns.R','EndTrialFuns.R','ExpFit.R','equipoiseFuns.R'), source)

args <- (commandArgs(TRUE)) ## load arguments from R CMD BATCH
print(args)
if(length(args)>0)  { ## Then cycle through each element of the list and evaluate the expressions.
    print(paste0('loading in ', args, ' from R CMD BATCH'))
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
    }
}else{
seed=120;trial="RCT";ord="TU";propInTrial=0.1;sdLogIndiv=1;delayUnit=7;immunoDelay=21;vaccEff=0.7;remStartFin="TRUE";remProtDel="TRUE";simNum=2880;batchdirnm="BigResults/gsRCT1";saveNm="gsRCT";nsims=1;reordLag=14;nboot=20;trialStartDate="2015-02-18"; gs=T
}

## hazType='Phenom'
## mu <- .15*yearToDays
## weeklyDecay <- .9
## cvWeeklyDecay <- .01
## cvClus <- .01
## cvClusTime <- .01
## gs <- T
## nsims <- 1#2000
## RCTendOption <- 3
## vaccEff <- 0
## doCFs <- T

seed=12;trial="RCT";gs="TRUE";ord="TU";propInTrial=0.05;delayUnit=7;immunoDelay=21;vaccEff=0.7;remStartFin="TRUE";remProtDel="TRUE";simNum=72;batchdirnm="BigResults/Equip1";saveNm="simSL-Equip-";nsims=170;reordLag=14;nboot=0;trialStartDate="2015-02-18"

verbose <- 1
parmArgs <- subsArgs(as.list(environment()), makeParms)
print(parmArgs)
parms <- do.call(makeParms, parmArgs)
saveFl <- file.path(batchdirnm, paste0(saveNm, sprintf("%06d", simNum),'.Rdata'))

system.time(sim <- simNtrials(seed=seed, parms=parms, N=nsims, verbFreq=10))
sim <- list(sim=sim, parms=parms, seed=seed, simNum=simNum)
save(sim, file = saveFl)

rm(list=ls(all=T))
gc()
