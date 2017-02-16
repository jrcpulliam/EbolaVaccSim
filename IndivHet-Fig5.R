if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevenbellan', Sys.info()['user'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('sbellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('nid', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/EbolaVaccSim/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T)); gc()
require(optiRum); require(RColorBrewer); require(boot); require(data.table); require(ggplot2); require(grid); require(reshape2); require(parallel)
sapply(c('multiplot.R','extractFXN.R','ggplotTheme.R'), source)
wid <- 10
heig <- 8
res <- 300
opacity <- 60

sapply(c('simFuns.R','AnalysisFuns.R','ExpFit.R'), source)

thing <- 'Equip-Fig5-v5'
figdir <- file.path('Figures', thing)
dir.create(figdir)
fls <- list.files('BigResults', pattern = paste0(thing,'-'), full.names=T)
colsv <- c("NT"='red', "SWCT"='orange', "VR"='dark green', "RCT-gs"='dodger blue', "RCT"='light blue', "RCT-rp"='purple', "RCT-gs-rp" = 'purple')

load(file.path('BigResults', paste0(thing, 'parmsMat','.Rdata')))
tvars <- c("trialStartDate", "propInTrial", "avHaz", "indivRRSeed" , "HazTrajSeed",'numClus','clusSuze','hazType','nbsize','mu','cvClus','cvClusTime','sdLogIndiv','weeklyDecay','cvWeeklyDecay','hazIntUnit')
tvars <- tvars[tvars %in% colnames(parmsMat)]
tpop <- unique(parmsMat[,tvars, with=F]) ## make sure to add other variables unique to population
tpop$tid <- 1:nrow(tpop)

irsk <- punq <- data.table()
for(ii in 1:length(fls)) {
    load(fls[ii])
    tmp <- merge(resList$punq,resList$SpopWorst,by = c('pid', 'propInTrial','lab'))
    tmp <- merge(tmp, resList$irsk[type=='marg', list(caseSpent = 6000*mean(spent_EV)), pid], by = 'pid')
    itmp <- resList$irsk
    itmp$tid <- ii ## **careful here**, need to make sure it always matches flow from parmsMat
    irsk <- rbind(irsk, itmp)
    speedpow <- with(resList, merge(parms[, list(pid, nbatch)], finTrials[,list(tcal, vaccEff,vaccGood, cvr,nbatch)], by = 'nbatch'))
    tmp <- merge(tmp, speedpow[, list(tcal = mean(tcal), .N), pid], by = 'pid')
    tmp$tid <- ii
    punq <- rbind(punq, tmp)
}
punq[,trialStartDate:=as.Date(trialStartDate)]
punq[,date:=format.Date(trialStartDate, '%b-%y')]
plunq <- punq[threshold==.05, list(lab, power, trialStartDate, threshold, above, above_EV, caseSpent, totCase, totCase_EV,avHaz, tcal,date)]
irsk <- merge(irsk, unique(punq[,.(tid,clusSize = as.numeric(clusSize))]), by='tid')## get clusSize in irsk for ordering purposes

if(!exists('hazT1')) {
    hazT1 <- setClusHaz(makePop(makeParms(trialStartDate="2014-10-01",propInTrial=0.05,avHaz="xTime",indivRRSeed=7,HazTrajSeed=7)))$hazT
    hazT1$avHaz <- 'xTime'
}
if(!exists('hazT2')) {
    hazT2 <- setClusHaz(makePop(makeParms(trialStartDate="2014-10-01",propInTrial=0.05,avHaz="",indivRRSeed=7,HazTrajSeed=7)))$hazT
    hazT2$avHaz <- ''
}
hazT <- rbind(hazT1, hazT2)

irsk[,ordX:=ordShowArm]
irsk[lab=='SWCT',ordX:=ordShow]

unique(parmsMat[,.(clusSize,avHaz, trialStartDate, .N), tid])
punq[threshold==.01,.(clusSize,avHaz, trialStartDate, lab,tid)]
names(irsk)

## Add a 6th of a cluster of space b/w each cluster for easier visualization
irsk[,ordX.sp:=ordX + floor((ordX-1)/clusSize)*(clusSize/6)]
avertableTab <- unique(irsk[lab %in% c('NT','VR'), .(arms = T, avertableRisk = inf[lab=='NT']-inf[lab=='VR'], ordX, ordX.sp, cVaccOrd, inf=inf[lab=='NT']), .(Oi,tid, clusSize)])
avertableTab2 <- merge(avertableTab[,.(Oi,tid,arms,avertableRisk,inf, clusSize, cVaccOrd)], unique(irsk[lab =='SWCT', .(Oi, ordX, ordX.sp)]), by ='Oi')
avertableTab2$arms <- F
avertableTab <- rbind(avertableTab, avertableTab2, fill=T)
                                        #avertableTab$ordShow <- NULL

avertableTab

####################################################################################################
## Plots
####################################################################################################
## Conditional on arms & order randomization
clusSizes <- irsk[,unique(clusSize)] 
clusSizes <- clusSizes[order(clusSizes)]

clusSizeDo <- clusSizes[1]

itmp <- irsk[clusSize==clusSizeDo]
tidDo <- itmp[,unique(tid)]
unique(itmp[,.(arm, armShown, type, lab, exmpl, tid)])
itmp <- itmp[tid==tidDo & ((arm==armShown & type=='cond' &  grepl('RCT',lab)) | (type=='condvd' & lab=='SWCT' & exmpl==T))]
itmp[, cols:=armShown]; itmp[cols=='cont',cols:='red']; itmp[cols=='vacc',cols:='dodger blue']
itmp <- itmp[lab!='RCT-rp']

atmp <- avertableTab[tid==tidDo & clusSize==clusSizeDo]

## for plotting
clusSizetmp <- punq[tid==tidDo, as.numeric(clusSize)][1]
Ntmp <- punq[tid==tidDo, as.numeric(clusSize)*as.numeric(numClus)][1]
mids <- seq(clusSizetmp/2, Ntmp-clusSizetmp/2, by = clusSizetmp)
mids <- mids + clusSizetmp/6*c(0:(length(mids)-1))
tcks <- rep(mids, each = 2) + c(-clusSizetmp/2,clusSizetmp/2)

numPanels <- 6
####################################################################################################
## INFECTION RISK, divided by randomization structure (for comparison to below)
pdf(file.path(figdir, paste0(clusSizeDo, '-', 'total infection risk by arm.pdf')), w = wid, h = heig)
xlim <- range(tcks)
par(mfrow=c(numPanels,1), mar = c(0,5,1,0), oma = c(4,1,0,0))
for(ll in itmp[,unique(lab)]) {
    itmp[lab==ll, plot(ordX.sp, inf, xlab='individual', ylab='', bty = 'n', type = 'h', col = cols, ylim = c(0,1), xlim=xlim, 
                       las = 1, xaxt='n', main = '')]
    if(ll == itmp[,unique(lab)][1])    legend('topright', col = c(itmp[,unique(cols)]), leg = c('control arm', 'vaccine arm'), pch = 15, cex = 1.3, bty = 'n')
    mtext(ll, side = 3, line = -2) 
}
mtext('cumulative infection risk without vaccination', 2, -1, outer = T)
axis(1, at = mids, lab = 1:20, lwd = 0)
dev.off()

####################################################################################################
## Averted by each trial with shadow of avertable, conditional on randomization assignment (exemplar for SWCT)
## CONDITIONAL
ylimavertable <- c(0, max(avertableTab[,avertableRisk]))
xlim <- c(0,Ntmp)
pdf(file.path(figdir, paste0(clusSizeDo, '-', 'avertable risk averted (conditional).pdf')), w = wid, h = heig)
par(mfrow=c(numPanels,1), mar = c(0,5,0,0), oma = c(5,2,0,0), ps = 15)
for(ll in itmp[,unique(lab)]) {
    atmp[arms==(ll!='SWCT'), plot(ordX.sp, avertableRisk, ylab ='', bty = 'n', type = 'h', col = 'gray', xlim = xlim, ylim = ylimavertable, xaxt='n', main = '', las = 1)]
    itmp[lab==ll, points(ordX.sp, avert_EV, xlab='individual', ylab='averted', bty = 'n', type = 'h', col = makeTransparent(cols,opacity))]
    if(ll==itmp[,unique(lab)][1]) {
        legend('topright', col = c('gray', itmp[,unique(cols)]), 
               leg = c('avertable risk', 'averted risk (control arm)', 'averted risk (vaccine arm)'), pch = 15, cex = 1.3, bty = 'n')
    }
    mtext(ll, 3, -4)
}
axis(1, at = mids, lab = 1:20, lwd = 0)
mtext(paste0('individuals by cluster (', clusSizetmp,' individuals per cluster)'), 1, outer=T, line = 3)
mtext('avertable risk', 2, ,outer=T, line = -1)
dev.off()

####################################################################################################
## Averted by each trial with shadow of avertable, Marginal on randomization assignment
## MARGINAL on arms & order randomization
itmpMarg <- irsk[tid==1 & type=='marg' & !lab %in% c('NT','VR')]
itmpMarg[, cols:=armShown]; itmpMarg[cols=='cont',cols:='red']; itmpMarg[cols=='vacc',cols:='dodger blue']
itmpMarg[lab=='SWCT', cols:='dodger blue']
itmpMarg <- itmpMarg[lab!='RCT-rp']

ylimavertable <- c(0, max(avertableTab[,avertableRisk]))
xlim <- c(0,Ntmp)
pdf(file.path(figdir, paste0(clusSizeDo, '-', 'avertable risk averted (marginal).pdf')), w = wid, h = heig)
par(mfrow=c(numPanels,1), mar = c(0,5,0,0), oma = c(5,2,0,0), ps = 15)
for(ll in itmpMarg[,unique(lab)]) {
    atmp[arms==(ll!='SWCT'), plot(ordX.sp, avertableRisk, ylab ='', bty = 'n', type = 'h', col = 'gray', xlim = xlim, ylim = ylimavertable, xaxt='n', main = '', las = 1)]
    itmpMarg[lab==ll, points(ordX.sp, avert_EV, xlab='individual', ylab='averted', bty = 'n', type = 'h', col = makeTransparent(cols,opacity))]
    if(ll==itmp[,unique(lab)][1]) legend('topright', col = c('gray', itmpMarg[,unique(cols)]), 
                                         leg = c('avertable risk', 'averted risk (control arm)', 'averted risk (vaccine arm)'), pch = 15, cex = 1.3, bty = 'n')
    mtext(ll, 3, -4)
}
axis(1, at = seq(0,Ntmp, by = clusSizetmp), lab = NA)
axis(1, at = seq(clusSizetmp/2,Ntmp-clusSizetmp/2, by = clusSizetmp), lab = 1:20, tick = 0)
mtext(paste0('individuals by cluster (', clusSizetmp,' individuals per cluster)'), 1, outer=T, line = 3)
mtext('avertable risk', 2, ,outer=T, line = -1)
dev.off()
####################################################################################################

####################################################################################################
## Selected clusters only
####################################################################################################
## Averted by each trial with shadow of avertable, conditional on randomization assignment (exemplar for SWCT)
## Zoom in on 1 cluster to explain plot

ylimavertable <- c(0, max(atmp[,avertableRisk]))
xlim <- c(0,Ntmp)

selClus <- 1#c(2,8,13)

itmp3 <- itmp[cVaccOrd %in% selClus]
itmp3 <- itmp3[order(ordX.sp)]
itmp3[,ord3:=1:length(ordX.sp), pid]
itmp3[,ord3.sp := ord3 + round(floor((ord3-1)/clusSizetmp)*clusSizetmp/6)]

atmp3 <- atmp[cVaccOrd %in% selClus]
atmp3 <- atmp3[order(ordX.sp)]
atmp3[,ord3:=1:length(ordX), arms]
atmp3[,ord3.sp := ord3 + round(floor((ord3-1)/clusSizetmp)*clusSizetmp/6)]

hazT$lwd <- 1
hazT[cluster %in% selClus, lwd:=5]
hazT[,col:=rainbow(length(unique(cluster)))[cluster]]
avHazDo <- punq[,unique(avHaz)]

pdf(file.path(figdir, paste0(clusSizeDo, '-', '3 clusters - avertable risk averted (conditional).pdf')), w = wid, h = heig)
every <- 1
xlim <- range(itmp3$ord3.sp)
par(mfrow=c(numPanels+1,1), mar = c(2,5,0,0), oma = c(5,2,0,0), ps = 16)
## highlight selected clusters
hazT[avHaz==avHazDo,plot(Date,clusHaz, lwd=0, bty = 'n', las = 1, ylab = 'hazard', yaxt = 'n', type = 'n')]
hazT[cluster!=selClus & avHaz==avHazDo, lines(Date, clusHaz, lwd = lwd, col = col), cluster]
hazT[cluster==selClus & avHaz==avHazDo, lines(Date, clusHaz, lwd = lwd, col = col), cluster]
#unique(hazT[, .(cluster,col)])[cluster %in% selClus, legend('topright', leg = cluster, col= col, bty = 'n', lwd = 3)]
par(mar=c(1,5,0,0))
## show risk breakdown
for(ll in itmp3[,unique(lab)]) {
    atmp3[arms==(ll!='SWCT') & ord3.sp%%every==0, 
          plot(ord3.sp, avertableRisk, ylab ='', bty = 'n', type = 'h', col = 'gray', xlim = xlim, ylim = ylimavertable, xaxt='n', main = '', las = 1, lwd = 2)]
    itmp3[lab==ll  & ord3.sp%%every==0, points(ord3.sp, avert_EV, xlab='individual', ylab='risk', bty = 'n', type = 'h', col = makeTransparent(cols,opacity), lwd=2)]
    if(ll==itmp[,unique(lab)][1]) legend('topright', col = c('gray', itmpMarg[,unique(cols)]), 
                                         leg = c('avertable', 'averted (control arm)', 'averted (vaccine arm)'), pch = 15, cex = 1.3, bty = 'n')
    mtext(ll, 3, -4)
}
axis(1, at = mids[1:length(selClus)], lab = paste0('cluster ', selClus), lwd = 0)
mtext('individuals', 1, outer=T, line = 3)
mtext('avertable risk', 2, ,outer=T, line = -1)
graphics.off()


