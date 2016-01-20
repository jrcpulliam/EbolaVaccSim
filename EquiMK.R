if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('wrangler', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/sbellan/wrangler/EbolaVaccSim/')
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

thing <- 'Equip-indivL'
batchdirnm <- file.path('BigResults',thing)
routdirnm <- file.path(batchdirnm,'Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
#tnms <- c('SWCT','RCT','FRCT')#,'CRCT')

numEach <- 24
nsims <- 85
print(paste0('doing ', nsims*numEach, ' per scenario with ', nsims , ' done on each of ', numEach, ' cores (batches)'))

## start dates
sdates <- seq.Date(as.Date('2014-10-01'), as.Date('2015-04-01'), by = 'month')
sdates <- sdates[1:length(sdates) %% 2 ==1]

## tnms <- c('RCT','SWCT','VR','NT')
ves <- NA
pits <- c(.025, .05, .1, .2)
parmsMatRCT <- as.data.table(expand.grid(
    batch =  1:numEach
  , trial = 'RCT'
  , gs = c(F,T)
  , ord = c('none','TU')
  , trialStartDate = sdates
  , propInTrial = pits
  , delayUnit = 7 ## c(0,7)
  , immunoDelay = c(21)
  , vaccEff = ves
  , randVaccProperties = T
  , vaccPropStrg='vaccProp1'
  , HazTrajSeed = 7
  , avHaz = c('', 'xTime','xClus','xClusxTime')
))
parmsMatSWCT <- as.data.table(expand.grid(
    batch =  1:numEach
  , trial = 'SWCT'
  , gs = F
  , ord = 'none'
  , trialStartDate = sdates
  , propInTrial = pits
  , delayUnit = 7 ## c(0,7)
  , immunoDelay = c(21)
  , vaccEff = ves
  , randVaccProperties = T
  , vaccPropStrg='vaccProp1'
  , HazTrajSeed = 7
  , avHaz = c('', 'xTime','xClus','xClusxTime')
))
parmsMatCFs <- as.data.table(expand.grid(
    batch =  1:numEach
  , trial = c('VR','NT')
  , gs = F
  , ord = 'TU'
  , trialStartDate = sdates
  , propInTrial = pits
  , delayUnit = 7 ## c(0,7)
  , immunoDelay = c(21)
  , vaccEff = ves
  , randVaccProperties = T
  , vaccPropStrg='vaccProp1'
  , HazTrajSeed = 7
  , avHaz = c('', 'xTime','xClus','xClusxTime')
))
parmsMat <- rbind(parmsMatRCT,parmsMatSWCT,parmsMatCFs)
parmsMat$returnEventTimes <- TRUE
parmsMat$StatsFxns <- 'doCoxME'
parmsMat[trial=='SWCT', StatsFxns:='doRelabel']
parmsMat[trial %in% c('NT','VR') , StatsFxns:=NA]
parmsMat <- parmsMat[order(avHaz,trial)]
parmsMat$rcmdbatch <- 1:nrow(parmsMat)
parmsMat$batchdirnm <- batchdirnm

nmtmp <- thing
parmsMat$saveNm <- nmtmp
parmsMat$nsims <- nsims
parmsMat$reordLag <- 14
parmsMat$nboot <- 200
nrow(parmsMat)

parmsMat[, simNumStart:=(batch-1)*nsims+1]
parmsMat[, simNumEnd:=(batch-1)*nsims+nsims]

parmsMat[order(gs), length(nboot), list(trial, ord, delayUnit, gs, vaccEff, avHaz, trialStartDate, propInTrial)]
nrow(parmsMat)
jbs <- NULL
jn <- 0

## batchdirnm <- file.path('BigResults',thing)
## fls <- list.files(batchdirnm, pattern=thing)
## fns <- as.numeric(sub('.Rdata','', sub(thing,'',fls)))
## parmsMatDo <- parmsMat[!rcmdbatch %in% fns]

## parmsMatDo <- parmsMat[trialStartDate %in% sdates[1:2] & trial=='RCT' & ord=='TU' & gs==T]
# parmsMatDo <- rbind(parmsMat[trial=='VR'][1],parmsMat[trial=='NT'][1],parmsMat[trial=='RCT'][1],parmsMat[trial=='SWCT'][1])
parmsMatDo <- parmsMat
nrow(parmsMatDo)
sink('SLsims.txt')
for(ii in parmsMatDo$rcmdbatch) {
    cmd <- "R CMD BATCH '--no-restore --no-save --args"
    cmd <- addParm(cmd, parmsMatDo, ii)
    cmd <- paste0(cmd, " ' startSim.R ", file.path(batchdirnm,'Routs', paste0(nmtmp, sprintf("%06d", ii),'.Rout')), 
                  sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
