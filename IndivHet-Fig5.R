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
    tmp <- merge(tmp, speedpow[, list(tcal = mean(tcal)), pid], by = 'pid')
    punq <- rbind(punq, tmp)
}
punq[,trialStartDate:=as.Date(trialStartDate)]
punq[,date:=format.Date(trialStartDate, '%b-%y')]
plunq <- punq[threshold==.05, list(lab, power, trialStartDate, threshold, above, above_EV, caseSpent, totCase, totCase_EV,avHaz, tcal,date)]

irsk[,ordX:=ordShowArm]
irsk[lab=='SWCT',ordX:=ordShow]

## add spacing b/w clusters for easier visualization
spc <- seq(0,50*19, by = 50)
spac <- data.table(ordX = 1:6000)
spac[, ordX.sp := ordX + rep(spc, each = 300)]
irsk <- merge(irsk, spac, by = 'ordX')

## for plotting
mids <- seq(150, 6000-150, by = 300)
mids <- mids + spc
tcks <- rep(mids, each = 2) + c(-150,150)

avertableTab <- unique(irsk[lab %in% c('NT','VR'), .(arms = T, avertableRisk = inf[lab=='NT']-inf[lab=='VR'], ordX, ordX.sp, cVaccOrd, inf=inf[lab=='NT']), .(Oi,tid)])
avertableTab2 <- merge(avertableTab[,.(Oi,tid,arms,avertableRisk,inf, cVaccOrd)], unique(irsk[lab =='SWCT', .(Oi, ordX, ordX.sp)]), by ='Oi')
avertableTab2$arms <- F
avertableTab <- rbind(avertableTab, avertableTab2, fill=T)
#avertableTab$ordShow <- NULL

####################################################################################################
## Plots
####################################################################################################
## Conditional on arms & order randomization
itmp <- irsk[tid==1 & ((arm==armShown & type=='cond' &  grepl('RCT',lab)) | (type=='condvd' & lab=='SWCT' & exmpl==T))]
itmp[, cols:=armShown]; itmp[cols=='cont',cols:='red']; itmp[cols=='vacc',cols:='dodger blue']
itmp <- itmp[lab!='RCT-rp']

####################################################################################################
## Look at individuals infection risk, but dividing them into their randomization structure (for comparison to below)
xlim <- range(tcks)
par(mfrow=c(6,1), mar = c(0,3,1,0), oma = c(4,1,0,0))
for(ll in itmp[,unique(lab)]) {
    itmp[lab==ll, plot(ordX.sp, inf, xlab='individual', ylab='risk spent', bty = 'n', type = 'h', col = cols, ylim = c(0,1), xlim=xlim, 
             las = 1, xaxt='n', main = '')]
    if(ll == itmp[,unique(lab)][1])    legend('topright', col = c(itmp[,unique(cols)]), leg = c('control arm', 'vaccine arm'), pch = 15, cex = 1.3, bty = 'n')
    mtext(ll, side = 3, line = -2) 
}
axis(1, at = mids, lab = 1:20, lwd = 0)

####################################################################################################
## Averted by each trial with shadow of avertable, conditional on randomization assignment (exemplar for SWCT)
ylimavertable <- c(0, max(avertableTab[,avertableRisk]))
xlim <- c(0,6000)
pdf(file.path(figdir, paste0('avertable risk averted (conditional).pdf')), w = wid, h = heig)
par(mfrow=c(6,1), mar = c(0,5,0,0), oma = c(5,2,0,0), ps = 15)
for(ll in itmp[,unique(lab)]) {
        avertableTab[arms==(ll!='SWCT'), plot(ordX.sp, avertableRisk, ylab ='', bty = 'n', type = 'h', col = 'gray', xlim = xlim, ylim = ylimavertable, xaxt='n', main = '', las = 1)]
        itmp[lab==ll, points(ordX.sp, avert_EV, xlab='individual', ylab='averted', bty = 'n', type = 'h', col = makeTransparent(cols,20))]
        if(ll==itmp[,unique(lab)][1]) {
            legend('topright', col = c('gray', itmp[,unique(cols)]), 
                   leg = c('avertable risk', 'averted risk (control arm)', 'averted risk (vaccine arm)'), pch = 15, cex = 1.3, bty = 'n')
        }
    mtext(ll, 3, -4)
}
axis(1, at = mids, lab = 1:20, lwd = 0)
mtext('individuals by cluster (300 individuals per cluster)', 1, outer=T, line = 3)
mtext('avertable risk', 2, ,outer=T, line = -1)
graphics.off()


####################################################################################################
## Averted by each trial with shadow of avertable, Marginal on randomization assignment

## Conditional on arms & order randomization
itmpMarg <- irsk[tid==1 & type=='marg' & !lab %in% c('NT','VR')]
itmpMarg[, cols:=armShown]; itmpMarg[cols=='cont',cols:='red']; itmpMarg[cols=='vacc',cols:='dodger blue']
itmpMarg[lab=='SWCT', cols:='dodger blue']
itmpMarg <- itmpMarg[lab!='RCT-rp']

ylimavertable <- c(0, max(avertableTab[,avertableRisk]))
xlim <- c(0,6000)
pdf(file.path(figdir, paste0('avertable risk averted (marginal).pdf')), w = wid, h = heig)
par(mfrow=c(6,1), mar = c(0,5,0,0), oma = c(5,2,0,0), ps = 15)
for(ll in itmpMarg[,unique(lab)]) {
        avertableTab[arms==(ll!='SWCT'), plot(ordX.sp, avertableRisk, ylab ='', bty = 'n', type = 'h', col = 'gray', xlim = xlim, ylim = ylimavertable, xaxt='n', main = '', las = 1)]
        itmpMarg[lab==ll, points(ordX.sp, avert_EV, xlab='individual', ylab='averted', bty = 'n', type = 'h', col = makeTransparent(cols,20))]
        if(ll==itmp[,unique(lab)][1]) legend('topright', col = c('gray', itmpMarg[,unique(cols)]), 
                   leg = c('avertable risk', 'averted risk (control arm)', 'averted risk (vaccine arm)'), pch = 15, cex = 1.3, bty = 'n')
    mtext(ll, 3, -4)
}
axis(1, at = seq(0,6000, by = 300), lab = NA)
axis(1, at = seq(150,6000-150, by = 300), lab = 1:20, tick = 0)
mtext('individuals by cluster (300 individuals per cluster)', 1, outer=T, line = 3)
mtext('avertable risk', 2, ,outer=T, line = -1)
graphics.off()
####################################################################################################

####################################################################################################
## 3 clusters only
####################################################################################################
## Averted by each trial with shadow of avertable, conditional on randomization assignment (exemplar for SWCT)
ylimavertable <- c(0, max(avertableTab[,avertableRisk]))
xlim <- c(0,6000)

selClus <- c(2,8,13)


itmp3 <- itmp[cVaccOrd %in% selClus]
itmp3 <- itmp3[order(ordX.sp)]
itmp3[,ord3:=1:length(ordX.sp), pid]
itmp3[,ord3.sp := ord3 + floor((ord3-1)/300)*50]


avertableTab3 <- avertableTab[cVaccOrd %in% selClus]
avertableTab3 <- avertableTab3[order(ordX.sp)]
avertableTab3[,ord3:=1:length(ordX), arms]
avertableTab3[,ord3.sp := ord3 + floor((ord3-1)/300)*50]

hazT$lwd <- 1
hazT[cluster %in% selClus, lwd:=5]
hazT[,col:=rainbow(length(unique(cluster)))[cluster]]

pdf(file.path(figdir, paste0('3 clusters - avertable risk averted (conditional).pdf')), w = wid, h = heig)
every <- 3
xlim <- range(itmp3$ord3.sp)
par(mfrow=c(7,1), mar = c(2,5,0,0), oma = c(5,2,0,0), ps = 16)
## highlight selected clusters
hazT[,plot(Date,clusHaz, lwd=0, bty = 'n', las = 1, ylab = 'hazard', yaxt = 'n', type = 'n')]
hazT[, lines(Date, clusHaz, lwd = lwd, col = col), cluster]
unique(hazT[, .(cluster,col)])[cluster %in% selClus, legend('topright', leg = cluster, col= col, bty = 'n', lwd = 3)]
par(mar=c(1,5,0,0))
## show risk breakdown
for(ll in itmp3[,unique(lab)]) {
    avertableTab3[arms==(ll!='SWCT') & ord3.sp%%every==1, 
                  plot(ord3.sp, avertableRisk, ylab ='', bty = 'n', type = 'h', col = 'gray', xlim = xlim, ylim = ylimavertable, xaxt='n', main = '', las = 1, lwd = 2)]
    itmp3[lab==ll  & ord3.sp%%every==1, points(ord3.sp, avert_EV, xlab='individual', ylab='averted', bty = 'n', type = 'h', col = makeTransparent(cols,150), lwd=2)]
    if(ll==itmp[,unique(lab)][1]) legend('topright', col = c('gray', itmpMarg[,unique(cols)]), 
           leg = c('avertable risk', 'averted risk (control arm)', 'averted risk (vaccine arm)'), pch = 15, cex = 1.3, bty = 'n')
    mtext(ll, 3, -4)
}
axis(1, at = mids, lab = paste0('cluster ', 1:20), lwd = 0)
mtext('individuals', 1, outer=T, line = 3)
mtext('avertable risk', 2, ,outer=T, line = -1)
graphics.off()


 
## ## Show total risk too (currently hard to display)
## xlim <- range(itmp3$ord3.sp)
## ylim <- c(0,1)
## lwd <- 2
## every <- 3
## par(mfrow=c(6,1), mar = c(0,5,0,0), oma = c(5,2,0,0), ps = 16)
## for(ll in itmp3[,unique(lab)]) {
##     avertableTab3[arms==(ll!='SWCT') & ord3.sp%%every==1, 
##                   plot(ord3.sp, inf, ylab ='', bty = 'n', type = 'h', col = gray(.3), xlim = xlim, ylim = ylim, xaxt='n', main = '', las = 1, lwd = lwd)]
##     avertableTab3[arms==(ll!='SWCT') & ord3.sp%%every==1, 
##                  points(ord3.sp, avertableRisk, col = gray(.8), lwd = lwd)]
##     itmp3[lab==ll  & ord3.sp%%every==1, points(ord3.sp, avert_EV, xlab='individual', ylab='averted', bty = 'n', type = 'h', col = makeTransparent(cols,150), lwd=lwd)]
##     if(ll==itmp[,unique(lab)][1]) legend('topright', col = c('gray', itmpMarg[,unique(cols)]), 
##            leg = c('avertable risk', 'averted risk (control arm)', 'averted risk (vaccine arm)'), pch = 15, cex = 1.3, bty = 'n')
##     mtext(ll, 3, -4)
## }
## axis(1, at = mids, lab = paste0('cluster ', 1:20), lwd = 0)
## mtext('individuals', 1, outer=T, line = 3)
## mtext('risk', 2, ,outer=T, line = -1)
## ##graphics.off()

## equipoise perturbed plot
## histogram within strata groups
## look at risk by person/strata by treatment assignment for trials
## vaccinating a greater % of people increases risk averted/spent, but still has problem of withholding treatment from individuals

## delay comparator doesn't provide much benefit when transmission is declining, or when risk accrual rate in the trial is sufficiently fast that control group will receive vaccine soon in a group sequential trial anyways. SO it really shows effect when risk is sustained/increasing, but not so high as to make the trial end very quickly, or if it's a very small trial
irsk[lab %in% c('RCT-gs-rp', 'RCT-gs-rp-cvd') & type == 'cond' & arm =='cont', diff(inf), Oi]

## show risk spent profiles in trial w/ placebo vs delayed vacc comparator, across range of trial sizes. Maybe show one very large trial population w/ lower risk, & one very small trial population w/ high risk but same amount of total hazard accumulation?
