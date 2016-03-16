if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('nid', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/EbolaVaccSim/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T)); gc()
require(optiRum); require(RColorBrewer); require(boot); require(data.table); require(ggplot2); require(grid); require(reshape2); require(parallel)

## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('extractFXN.R'), source)

thing <- 'Equip-RRcat'
## Load VaccProp & hazT
load('data/vaccProp1.Rdata')
vaccProp <- vaccProp1
vaccProp[, simNum:=1:length(vaccEff)]
## Make Figure Folder
figdir <- file.path('Figures', thing)
dir.create(figdir)

## batchdirnm <- file.path('BigResults',thing)
## fls <- list.files(batchdirnm, pattern=thing, full.names = T)
## ## fls <- fls[grepl(2305,fls)]
## eos <- extractOneSim(fls[1], indivLev=T, verbose = 0)

load(file.path('BigResults', paste0(thing, 'parmsMat','.Rdata')))
## start dates
sdates <- seq.Date(as.Date('2014-10-01'), as.Date('2015-04-01'), by = 'month')
sdates <- sdates[1:length(sdates) %% 2 ==1]

unique(parmsMat[avHaz=='xTime' & propInTrial==c(.05) & trialStartDate==c('2014-10-01'), list(avHaz, tid)])
tidsDo <- unique(parmsMat[propInTrial == c(.05) & trialStartDate %in% sdates[c(1,3)] & avHaz %in% c('', 'xTime'), tid] )

for(ti in tidsDo) {
    print(ti)
    procAll(tidDo = ti,0)
}

for(ti in tidsDo) {
    load(file=file.path('BigResults',paste0(thing, '-', tidDo, '.Rdata')))
}
