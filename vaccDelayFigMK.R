if(grepl('sbellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebe', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('nid', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/EbolaVaccSim/')
if(grepl('wrangler', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/sbellan/wrangler/EbolaVaccSim/')
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R'), source)
## CHANGE SWCT to relabel when want power again ***

thing <- 'Equip-Fig-SX-vaccDel-3cat'
batchdirnm <- file.path('BigResults',thing)
routdirnm <- file.path(batchdirnm,'Routs')
if(!file.exists(batchdirnm)) dir.create(batchdirnm)
if(!file.exists(routdirnm)) dir.create(routdirnm)
#tnms <- c('SWCT','RCT','FRCT')#,'CRCT')

numEach <- 24
nsims <- 85
print(paste0('doing ', nsims*numEach, ' per scenario with ', nsims , ' done on each of ', numEach, ' cores (batches)'))

## start dates
sdates <- as.Date('2014-10-01')

ves <- NA
pits <- .05
avHazs <- c('', 'xTime')
clusSizes <- 40*c(1:5)
parmsMatRCT <- as.data.table(expand.grid(
    batch =  1:numEach
  , clusSize = clusSizes
  , nsims = nsims
  , trial = 'RCT'
  , gs = T
  , ord = 'TU'
  , contVaccDelay = c(NA, 7*9)
  , maxRRcat = c(0, 25)
  , trialStartDate = sdates
  , propInTrial = pits
  , vaccEff = ves
  , avHaz = avHazs
))
## only do the threshold and vaccContDelay for gsTU trials
parmsMatRCT <- parmsMatRCT[!(!is.na(contVaccDelay) & (gs==F | ord=='none'))]
parmsMatRCT <- parmsMatRCT[!(maxRRcat>0 & (gs==F | ord=='none'))]
parmsMatRCT <- parmsMatRCT[!(maxRRcat>0 & !is.na(contVaccDelay))]
parmsMatCFs <- as.data.table(expand.grid(
    batch =  1:numEach ## do for average vacc eff
  , clusSize = clusSizes
  , nsims = nsims ## 
  , trial = c('VR','NT')
  , gs = F
  , contVaccDelay = NA
  , maxRRcat = 0
  , ord = 'TU'
  , trialStartDate = sdates
  , propInTrial = pits
  , vaccEff = ves
  , avHaz = avHazs
))

parmsMat <- rbind(parmsMatRCT,parmsMatCFs)
parmsMat <- within(parmsMat, { vaccPropStrg='vaccProp1'; HazTrajSeed=7; indivRRSeed=7; returnEventTimes=TRUE; immunoDelay=21; delayUnit=7; randVaccProperties=T;
                               DoIndivRRcat=T})
parmsMat$StatsFxns <- 'doCoxME'
parmsMat[trial %in% c('NT','VR') , StatsFxns:=NA]
parmsMat <- parmsMat[order(avHaz,trial)]
parmsMat$rcmdbatch <- 1:nrow(parmsMat)
parmsMat$batchdirnm <- batchdirnm
parmsMat$numClus <- 10 #### careful, change if necessary!
parmsMat$constRiskXClusSize <- T
names(parmsMat)

## variables that specify a trial population
tvars <- c("trialStartDate", "propInTrial", "avHaz", "indivRRSeed" , "HazTrajSeed",'numClus','clusSize','hazType','nbsize','mu','cvClus','cvClusTime','sdLogIndiv','weeklyDecay','cvWeeklyDecay','hazIntUnit')
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
unique(parmsMat[,.(clusSize,avHaz, trialStartDate, .N), tid])

nmtmp <- thing
parmsMat <- within(parmsMat, {saveNm = nmtmp; reordLag = 14; nboot = 200})
nrow(parmsMat)

parmsMat[, simNumStart:=(batch-1)*nsims+1]
parmsMat[, simNumEnd:=(batch-1)*nsims+nsims]
save(parmsMat, file=file.path('BigResults', paste0(thing, 'parmsMat','.Rdata')))

## tidsDo <- tpop[propInTrial == c(.05) & trialStartDate %in% sdates[c(1,3)] & avHaz %in% c('', 'xTime'), tid]
## tidsDo <- tpop[propInTrial == c(.05) & avHaz %in% c('', 'xTime'), tid]
##tidsDo <- tpop[propInTrial %in% c(.05,.1) & trialStartDate %in% c('2014-10-01','2014-12-01'), tid]
## tidsDo <- tpop[propInTrial == c(.05) & avHaz %in% '' & trialStartDate=='2014-10-01', tid] ### CHANGE*** as appropriate
## parmsMatDo <- parmsMat[tid %in% tidsDo]

parmsMatDo <- parmsMat##[propInTrial == c(.05) & avHaz %in% '' & trialStartDate=='2014-10-01']
nrow(parmsMatDo)

## fls <- list.files(file.path('BigResults', thing), pattern='.Rdata')
## rcmdsDone <- as.numeric(gsub('.Rdata', '', gsub(thing, '', fls)))
## rcmdsDone <- rcmdsDone[order(rcmdsDone)]
## length(rcmdsDone)
## toDo <- parmsMat[!(rcmdbatch %in% rcmdsDone), rcmdbatch]
## length(toDo)
## parmsMat[rcmdbatch %in% toDo, table(trial)]
## parmsMat[!rcmdbatch %in% toDo, table(trial)]
## parmsMat[,table(trial)]
## parmsMatDo <- parmsMat[rcmdbatch %in% toDo]

parmsMatDo[,table(trial)]

## summary of runs
unique(parmsMatDo[,.(nruns = nsims *.N),.(clusSize,avHaz, trialStartDate, gs, ord, maxRRcat, contVaccDelay, trial, tid)])

#parmsMatDo <- parmsMatDo[trial %in% c('NT','VR')]

##parmsMatDo <- parmsMatDo[!trial %in% c('NT','VR','SWCT')]
nrow(parmsMatDo)
sink('SLsimsCVD.txt')
for(ii in parmsMatDo$rcmdbatch) {
    cmd <- "R CMD BATCH '--no-restore --no-save --args"
    cmd <- addParm(cmd, parmsMatDo, ii)
    cmd <- paste0(cmd, " ' startSim.R ", file.path(batchdirnm,'Routs', paste0(nmtmp, sprintf("%06d", ii),'.Rout')), 
                  sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
