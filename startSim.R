if(grepl("Stevens-MacBook-Pro.local", Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('nid', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/EbolaVaccSim/')
if(grepl('wrangler', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/sbellan/wrangler/EbolaVaccSim/')
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
## pid=38;tid=21;batch=1;rcmdbatch=2385;trialStartDate="2014-10-01";propInTrial=0.05;avHaz="xTime";indivRRSeed=7;HazTrajSeed=7;nsims=5;trial="RCT";gs="FALSE";ord="TU";vaccEff="NA";randVaccProperties="TRUE";delayUnit=7;immunoDelay=21;returnEventTimes="TRUE";vaccPropStrg="vaccProp1";StatsFxns="doCoxME";batchdirnm="BigResults/Equip-irskHybrid";nboot=200;reordLag=14;saveNm="Equip-irskHybrid";simNumStart=1;simNumEnd=85; DoIndivRRcat=T

pid=6;tid=1;batch=1;rcmdbatch=121;trialStartDate="2014-12-01";propInTrial=0.05;avHaz="";indivRRSeed=7;HazTrajSeed=7;nsims=2;trial="RCT";gs="TRUE";ord="TU";contVaccDelay=7*8;maxRRcat=25;vaccEff=.8;DoIndivRRcat="TRUE";randVaccProperties="FALSE";delayUnit=7;immunoDelay=21;returnEventTimes="TRUE";vaccPropStrg="vaccProp1";StatsFxns="doCoxME";batchdirnm="BigResults/Equip-Fig5-v1";nboot=200;reordLag=14;saveNm="Equip-Fig5-v1";simNumStart=1;simNumEnd=2

pid=9;tid=1;batch=24;rcmdbatch=408;trialStartDate="2014-12-01";propInTrial=0.05;avHaz="";indivRRSeed=7;HazTrajSeed=7;nsims=2;trial="NT";gs="FALSE";ord="TU";contVaccDelay=NA;maxRRcat=0;vaccEff="NA";DoIndivRRcat="TRUE";randVaccProperties="TRUE";delayUnit=7;immunoDelay=21;returnEventTimes="TRUE";vaccPropStrg="vaccProp1";StatsFxns="NA";batchdirnm="BigResults/Equip-Fig5-v4";nboot=200;reordLag=14;saveNm="Equip-Fig5-v4";simNumStart=1956;simNumEnd=1958

## pid=7;tid=1;batch=1;rcmdbatch=289;trialStartDate="2014-12-01";propInTrial=0.05;avHaz="";indivRRSeed=7;HazTrajSeed=7;nsims=1;trial="SWCT";gs="FALSE";ord="none";contVaccDelay=7*8;maxRRcat=30;vaccEff="NA";DoIndivRRcat="TRUE";randVaccProperties="TRUE";delayUnit=7;immunoDelay=21;returnEventTimes="TRUE";vaccPropStrg="vaccProp1";StatsFxns="doRelabel";batchdirnm="BigResults/Equip-Fig5-v2";nboot=200;reordLag=14;saveNm="Equip-Fig5-v2";simNumStart=1;simNumEnd=1

}
load('data/vaccProp1.Rdata')

## vaccPropStrg <- 'vaccProp1'
if(randVaccProperties)
   vaccProp <- get(vaccPropStrg) else vaccProp <- NA

verbose <- 1
parmArgs <- subsArgs(as.list(environment()), makeParms)
print(parmArgs)
parms <- do.call(makeParms, parmArgs)
saveFl <- file.path(batchdirnm, paste0(saveNm, sprintf("%06d", rcmdbatch),'.Rdata'))
simArgs <- list(batch=batch, parms=parms, N=nsims, verbFreq=10, vaccProp=vaccProp, returnEventTimes=returnEventTimes,
                simNums=simNumStart:simNumEnd)

hazT <- NA
if(!is.na(HazTrajSeed)) {
#    rnorm(1)
    hazT <- setClusHaz(makePop(parms))$hazT
    ## pdf('Figures/tst.pdf')
    ## ggplot(hazT) + geom_line(aes(Date,clusHaz, col=factor(cluster)))
    ## dev.off()
    ##    save(hazT, file=paste0('BigResults/Equip-indivL/hazT',HazTrajSeed,'.Rdata'))
}

system.time(sim <- do.call(simNtrials, simArgs))
sim <- list(sim=sim, parms=parms, batch=batch, rcmdbatch=rcmdbatch, hazT=hazT)
save(sim, file = saveFl)

rm(list=ls(all=T))
gc()

## * analytics for CFR & equipoise with cfs, just use daly multiplier on sae with alread run simulations (cfr increases divide and reduces overlap with perfect equipoise. inf risk does same but also affects power)

## * size of trial decreases each individuals' anticipated equipoise dilemma

## * herd immunity (sherif)

## * ring vaccination

## * CRCT

## JP notes
## provide people a tool to draw a curve infection spending tool themselves and then compare different designs. state assumptions & turn it into a curve (shiny)

## LAM notes

## calculate power / inf spent for all trials, plot that over scenarios; might need to show absolute power too still
## excess risk taken on distribution, ppl above a certain threshold should be given a choice regardless of info provided
## informed consent, choice to be given experimental vaccine vs being randomized, working in altruism
## what does a trialist do with our framework, how can they decide who is at what level? simple BMJ version with suspected inf risk & CFR
