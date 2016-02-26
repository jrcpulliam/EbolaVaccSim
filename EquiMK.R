if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('nid', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/EbolaVaccSim/')
if(grepl('wrangler', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/sbellan/wrangler/EbolaVaccSim/')
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)
## CHANGE SWCT to relabel when want power again ***

thing <- 'Equip-irskHybrid'
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
    , nsims = nsims
  , trial = 'RCT'
  , gs = c(F,T)
  , ord = c('none','TU')
  , trialStartDate = sdates
  , propInTrial = pits
  , vaccEff = ves
  , avHaz = c('', 'xTime','xClus','xClusxTime')
))
parmsMatSWCT <- as.data.table(expand.grid(
    batch =  1:numEach
  , nsims = nsims
  , trial = 'SWCT'
  , gs = F
  , ord = 'none'
  , trialStartDate = sdates
  , propInTrial = pits
  , vaccEff = ves
  , avHaz = c('', 'xTime','xClus','xClusxTime')
))
parmsMatCFs <- as.data.table(expand.grid(
    batch =  1 ## only one for CFs
  , nsims = 1 ## 
  , trial = c('VR','NT')
  , gs = F
  , ord = 'TU'
  , trialStartDate = sdates
  , propInTrial = pits
  , vaccEff = ves
  , avHaz = c('', 'xTime','xClus','xClusxTime')
))

parmsMat <- rbind(parmsMatRCT,parmsMatSWCT,parmsMatCFs)
parmsMat <- within(parmsMat, { vaccPropStrg='vaccProp1'; HazTrajSeed=7; indivRRSeed=7; returnEventTimes=TRUE; immunoDelay=21; delayUnit=7; randVaccProperties=T})
parmsMat$StatsFxns <- 'doCoxME'
parmsMat[trial=='SWCT', StatsFxns:='doCoxME'] ### CHANGE BACK
parmsMat[trial %in% c('NT','VR') , StatsFxns:=NA]
parmsMat <- parmsMat[order(avHaz,trial)]
parmsMat$rcmdbatch <- 1:nrow(parmsMat)
parmsMat$batchdirnm <- batchdirnm

names(parmsMat)

## variables that specify a trial population
tvars <- c("trialStartDate", "propInTrial", "avHaz", "indivRRSeed" , "HazTrajSeed",'numClus','clusSuze','hazType','nbsize','mu','cvClus','cvClusTime','sdLogIndiv','weeklyDecay','cvWeeklyDecay','hazIntUnit')
tvars <- tvars[tvars %in% colnames(parmsMat)]
tpop <- unique(parmsMat[,tvars, with=F]) ## make sure to add other variables unique to population
tpop <- data.table(tid=1:nrow(tpop), tpop)
parmsMat <- merge(tpop, parmsMat, by=names(tpop)[names(tpop)!='tid'])
frontcols <- c('tid','batch')
setcolorder(parmsMat, c(frontcols, names(parmsMat)[!names(parmsMat) %in% frontcols]))
## variables that specify a parameter combination (which may be spread over batches across cores)
punq <- parmsMat[,!c('tid','batch','rcmdbatch'),with=F]
setkey(punq, NULL)
punq <- unique(punq)
punq <- data.table(pid = 1:nrow(punq), punq)
setkey(punq, pid)
parmsMat <- merge(punq, parmsMat, by = names(punq)[names(punq)!='pid'])
setkey(parmsMat, tid, pid, batch)
frontcols <- c('pid','tid','batch', 'rcmdbatch')
setcolorder(parmsMat, c(frontcols, names(parmsMat)[!names(parmsMat) %in% frontcols]))
parmsMat

nmtmp <- thing
parmsMat <- within(parmsMat, {saveNm = nmtmp; reordLag = 14; nboot = 200})
nrow(parmsMat)

parmsMat[, simNumStart:=(batch-1)*nsims+1]
parmsMat[, simNumEnd:=(batch-1)*nsims+nsims]
save(parmsMat, file=file.path('BigResults', paste0(thing, 'parmsMat','.Rdata')))

tidsDo <- tpop[propInTrial == c(.05) & trialStartDate == c('2014-10-01') & avHaz =='xTime', tid]
##tidsDo <- tpop[propInTrial %in% c(.05,.1) & trialStartDate %in% c('2014-10-01','2014-12-01'), tid]

parmsMatDo <- parmsMat[tid %in% tidsDo]

##parmsMatDo <- parmsMatDo[!trial %in% c('NT','VR','SWCT')]
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

## batchdirnm <- file.path('BigResults',thing)
## fls <- list.files(batchdirnm, pattern=thing)
## fns <- as.numeric(sub('.Rdata','', sub(thing,'',fls)))
## parmsMatDo <- parmsMat[!rcmdbatch %in% fns] ## 174
