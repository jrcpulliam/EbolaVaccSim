if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('wrangler', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/sbellan/wrangler/EbolaVaccSim/')
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)

doCFs <- F

thing <- paste0('Equip-ByTrialDate', 'CFs'[doCFs])
batchdirnm <- file.path('BigResults',thing)
routdirnm <- file.path(batchdirnm,'Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
#tnms <- c('SWCT','RCT','FRCT')#,'CRCT')
tnms <- c('RCT','SWCT')
numEach <- 24
nsims <- 85
print(paste0('doing ', nsims*numEach, ' per scenario with ', nsims , ' done on each of ', numEach, ' cores (batches)'))

## start dates
sdates <- seq.Date(as.Date('2014-10-01'), as.Date('2015-04-01'), by = 'month')
sdates <- sdates[1:length(sdates) %% 2 ==1]
## sdates <- '2015-02-18'


if(doCFs) { 
    gs <- F
    ord <- 'TU'
    tnms <- 'RCT'
}else{
    gs <- c(F, T)
    ord <- c('none','TU')
}
ves <- NA
pits <- c(.025, .05, .1, .2)
parmsMat <- as.data.table(expand.grid(
    batch =  1:numEach
  , trial = tnms
  , gs = gs
  , ord = ord
  , trialStartDate = sdates
  , propInTrial = pits
  , delayUnit = 7 ## c(0,7)
  , immunoDelay = c(21)
  , vaccEff = ves
  , randVaccProperties = T
  , vaccPropStrg='vaccProp1'
  , numCFs = 1
  , HazTrajSeed = 7
    , avHaz = c('', 'xTime','xClus','xClusxTime')
))
parmsMat$returnEventTimes <- TRUE
parmsMat$doCFs <- doCFs
parmsMat <- parmsMat[order(trial)]
parmsMat$StatsFxns <- 'doCoxME'
parmsMat[trial=='SWCT', StatsFxns:='doRelabel']
parmsMat <- parmsMat[!(trial=='SWCT' & (delayUnit==0 | ord=='TU'))] ## SWCT must have delay and cannot be ordered
parmsMat <- parmsMat[!(trial=='SWCT' & gs)] ## SWCT must have delay and cannot be ordered
parmsMat <- parmsMat[!(delayUnit==0 & ord=='TU')] ## ordering is meaningless with simultaneous instant vacc
parmsMat <- parmsMat[ !(delayUnit==0 & trial=='FRCT')]  ## FRCT = RCT when delayUnit=0
parmsMat <- parmsMat[order(avHaz)]

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


batchdirnm <- file.path('BigResults',thing)
fls <- list.files(batchdirnm, pattern=thing)
fns <- as.numeric(sub('.Rdata','', sub(thing,'',fls)))

parmsMatDo <- parmsMat[!rcmdbatch %in% fns]
nrow(parmsMatDo)
sink(paste0('SLsims','CFs'[doCFs],'.txt'))
for(ii in parmsMatDo$rcmdbatch) {
    cmd <- "R CMD BATCH '--no-restore --no-save --args"
    cmd <- addParm(cmd, parmsMatDo, ii)
    cmd <- paste0(cmd, " ' startSim.R ", file.path(batchdirnm,'Routs', paste0(nmtmp, sprintf("%06d", ii),'.Rout')), 
                  sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
