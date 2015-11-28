if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('wrangler', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/sbellan/wrangler/EbolaVaccSim/')
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)
 
thing <- 'Equip-randCFs'
batchdirnm <- file.path('BigResults',thing)
routdirnm <- file.path(batchdirnm,'Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
tnms <- c('SWCT','RCT','FRCT')#,'CRCT')
tnms <- 'RCT'
numEach <- 40

ves <- c(0,.7)
pits <- c(.05)
parmsMat <- as.data.table(expand.grid(
    batch =  1:numEach
    , trial = tnms
  , gs = c(F, T)
  , ord = c('none','TU')
  , propInTrial = pits
  , delayUnit = c(0,7)
  , immunoDelay = c(21)
  , vaccEff = ves
  , randVaccProperties = T
  , vaccPropStrg='vaccProp1'
  , numCFs = 500
))
parmsMat$remStartFin <- TRUE ##***
parmsMat$remProtDel <- TRUE
parmsMat$returnEventTimes <- TRUE
parmsMat[,doCFs:= numCFs>0]
parmsMat <- parmsMat[!(trial=='SWCT' & (delayUnit==0 | ord=='TU'))] ## SWCT must have delay and cannot be ordered
parmsMat <- parmsMat[!(trial=='SWCT' & gs)] ## SWCT must have delay and cannot be ordered
parmsMat <- parmsMat[!(delayUnit==0 & ord=='TU')] ## ordering is meaningless with simultaneous instant vacc
parmsMat <- parmsMat[ !(delayUnit==0 & trial=='FRCT')]  ## FRCT = RCT when delayUnit=0
parmsMat$rcmdbatch <- 1:nrow(parmsMat)
parmsMat$batchdirnm <- batchdirnm
nmtmp <- thing
parmsMat$saveNm <- nmtmp
parmsMat$nsims <- 25 ## 85*24 is ~ 2000 simulations each (2040 but we'll round)
parmsMat$reordLag <- 14
parmsMat$nboot <- 0
parmsMat$trialStartDate <- '2015-02-18'
nrow(parmsMat)

parmsMat[, simNumStart:=(batch-1)*nsims+1]
parmsMat[, simNumEnd:=(batch-1)*nsims+nsims]

parmsMat[order(gs), length(nboot), list(trial, ord, delayUnit, gs, vaccEff)]
nrow(parmsMat)
jbs <- NULL
jn <- 0

parmsMatDo <- parmsMat
sink(paste0('SLsims.txt'))
for(ii in parmsMatDo$simNum) {
    cmd <- "R CMD BATCH '--no-restore --no-save --args"
    cmd <- addParm(cmd, parmsMatDo, ii)
    cmd <- paste0(cmd, " ' startSim.R ", file.path(batchdirnm,'Routs', paste0(nmtmp, sprintf("%06d", ii),'.Rout')), 
                  sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
