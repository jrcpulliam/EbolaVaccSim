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
batch=1;trial="RCT";ord="TU";propInTrial=0.1;sdLogIndiv=1;delayUnit=7;immunoDelay=21;vaccEff=0.7;remStartFin="TRUE";remProtDel="TRUE";rcmdbatch=2880;batchdirnm="BigResults/gsRCT1";saveNm="gsRCT";nsims=1;reordLag=14;nboot=20;trialStartDate="2015-02-18"; gs=T;randVaccProperties=F;vaccPropStrg=NA;simNumStart=(batch-1)*nsims+1;simNumEnd=(batch-1)*nsims+nsims;returnEventTimes=T;pSAE=10^-2;HazTrajSeed=1
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

system.time(sim <- simNtrialsWRP(simArgs))
sim <- list(sim=sim, parms=parms, batch=batch, rcmdbatch=rcmdbatch)
save(sim, file = saveFl)

rm(list=ls(all=T))
gc()
