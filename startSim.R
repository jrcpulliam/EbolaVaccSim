if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
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
batch=1;trial="SWCT";gs="FALSE";ord="none";trialStartDate="2015-04-01";propInTrial=0.05;delayUnit=7;immunoDelay=21;vaccEff="NA";randVaccProperties="TRUE";vaccPropStrg="vaccProp1";HazTrajSeed=7;avHaz="xTime";returnEventTimes="TRUE";StatsFxns="doCoxME";rcmdbatch=4393;batchdirnm="BigResults/Equip-indivL";saveNm="Equip-indivL";nsims=3;reordLag=14;nboot=200;simNumStart=1;simNumEnd=85;indivRRSeed=1;
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

## * stratify equip by risk (must returnEventTimes for all sims, not just cf; or some aggregate version)

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
