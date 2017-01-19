if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
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

thing <- 'Equip-Fig5-v4'
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

## add spacing b/w clusters for easier visualization
spc <- seq(0,50*19, by = 50)
spac <- data.table(ordShowArm = 1:6000)
spac[, ordShowArm.sp := ordShowArm + rep(spc, each = 300)]
irsk <- merge(irsk, spac, by = 'ordShowArm')
spac2 <- data.table(ordShow = 1:6000)
spac2[, ordShow.sp := ordShow + rep(spc, each = 300)]
irsk <- merge(irsk, spac2, by = 'ordShow')

## for plotting
mids <- seq(150, 6000-150, by = 300)
mids <- mids + spc
tcks <- rep(mids, each = 2) + c(-150,150)

avertableTab <- irsk[lab %in% c('NT','VR'), .(avertableRisk = inf[lab=='NT']-inf[lab=='VR'], ordShow, ordShowArm), .(Oi,tid)]
avertableTab <- merge(avertableTab, spac, by = 'ordShowArm')
avertableTab <- merge(avertableTab, spac2, by = 'ordShow')

####################################################################################################
## Plots
####################################################################################################
## Conditional on arms & order randomization
itmp <- irsk[tid==1 & ((arm==armShown & type=='cond' &  grepl('RCT',lab)) | (type=='condvd' & lab=='SWCT' & exmpl==T))]
itmp[, cols:=armShown]; itmp[cols=='cont',cols:='red']; itmp[cols=='vacc',cols:='dodger blue']
itmp <- itmp[lab!='RCT-rp']

####################################################################################################
## Look at individuals infection risk, but dividing them into their randomization structure (for comparison to below)
xlim <- c(0,6000)
par(mfrow=c(6,1), mar = c(0,3,1,0), oma = c(4,1,0,0))
for(ll in itmp[,unique(lab)]) {
    itmp[lab==ll, plot(ordShowArm.sp, inf, xlab='individual', ylab='risk spent', bty = 'n', type = 'h', col = cols, ylim = c(0,1), xlim=xlim, 
             las = 1, xaxt='n', main = '')]
    mtext(ll, side = 3, line = -2) 
}
axis(1, at = mids, lab = 1:20, lwd = 0)

####################################################################################################
## Averted by each trial with shadow of avertable, conditional on randomization assignment (exemplar for SWCT)
ylimavertable <- c(0, max(avertableTab[,avertableRisk]))
xlim <- c(0,6000)
pdf(file.path(figdir, paste0('avertable risk averted (conditional).pdf')), w = wid, h = heig)
step <- 0
par(mfrow=c(6,1), mar = c(0,5,0,0), oma = c(5,2,0,0), ps = 15)
for(ll in itmp[,unique(lab)]) {
    if(ll!='SWCT'){
        avertableTab[, plot(ordShowArm.sp, avertableRisk, ylab ='', bty = 'n', type = 'h', col = 'gray', xlim = xlim, ylim = ylimavertable, xaxt='n', main = '', las = 1)]
        itmp[lab==ll, points(ordShowArm.sp, avert_EV, xlab='individual', ylab='averted', bty = 'n', type = 'h', col = makeTransparent(cols,20))]
        if(step==0){
            legend('topright', col = c('gray', itmp[,unique(cols)]), 
                   leg = c('avertable risk', 'averted risk (control arm)', 'averted risk (vaccine arm)'), pch = 15, cex = 1.3, bty = 'n')
            step <- 1
        }
    }else{
        avertableTab[, plot(ordShow.sp, avertableRisk, ylab ='', bty = 'n', type = 'h', col = 'gray', xlim = xlim, ylim = ylimavertable, xaxt='n', main = '', las = 1)]
        itmp[lab==ll, points(ordShow.sp, avert_EV, xlab='individual', ylab='averted', bty = 'n', type = 'h', col = makeTransparent(cols,20))]
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
step <- 0
pdf(file.path(figdir, paste0('avertable risk averted (marginal).pdf')), w = wid, h = heig)
par(mfrow=c(6,1), mar = c(0,5,0,0), oma = c(5,2,0,0), ps = 15)
for(ll in itmpMarg[,unique(lab)]) {
    if(ll!='SWCT'){
        avertableTab[, plot(ordShowArm.sp, avertableRisk, ylab ='', bty = 'n', type = 'h', col = 'gray', xlim = xlim, ylim = ylimavertable, xaxt='n', main = '', las = 1)]
        itmpMarg[lab==ll, points(ordShowArm.sp, avert_EV, xlab='individual', ylab='averted', bty = 'n', type = 'h', col = makeTransparent(cols,20))]
        if(step==0){
            legend('topright', col = c('gray', itmpMarg[,unique(cols)]), 
                   leg = c('avertable risk', 'averted risk (control arm)', 'averted risk (vaccine arm)'), pch = 15, cex = 1.3, bty = 'n')
            step <- 1
        }
    }else{
        avertableTab[, plot(ordShow.sp, avertableRisk, ylab ='', bty = 'n', type = 'h', col = 'gray', xlim = xlim, ylim = ylimavertable, xaxt='n', main = '', las = 1)]
        itmpMarg[lab==ll, points(ordShow.sp, avert_EV, xlab='individual', ylab='averted', bty = 'n', type = 'h', col = makeTransparent(cols,20))]
    }
    mtext(ll, 3, -4)
}
axis(1, at = seq(0,6000, by = 300), lab = NA)
axis(1, at = seq(150,6000-150, by = 300), lab = 1:20, tick = 0)
mtext('individuals by cluster (300 individuals per cluster)', 1, outer=T, line = 3)
mtext('avertable risk', 2, ,outer=T, line = -1)
graphics.off()
####################################################################################################

par(mfrow=c(6,1), mar = c(0,3,1,0), oma = c(1,1,0,0))
for(ll in mtmp[,unique(lab)]) {
    mtmp[lab==ll, plot(ord, spent_EV, xlab='individual', ylab='risk spent', bty = 'n', type = 'h', col = cols, 
             las = 1, xaxt='n', main =ll)]
    with(mtmp[lab==ll], points(ord, -avert_EV, type = 'h', col = makeTransparent(cols, alpha = 250)))
    abline(h=0, lty = 1)
}
##graphics.off()

## look at all yvars for each trial type to see if things look right
pdf(file.path(figdir, paste0('irsk InfSpentAvert by trial.pdf')), w = wid, h = heig)
for(ll in mtmp[,unique(lab)]) {
    par(mfrow=c(6,1), mar = c(0,3,1,0), oma = c(1,1,0,0))
    for(yvar in c('inf', 'avert','avert_EV', 'spent', 'spent_EV')) {
        ymax <- ifelse(yvar=='inf', 1, .15)
        mtmp[lab==ll, plot(ord, get(yvar), xlab='individual', ylab='risk spent', bty = 'n', type = 'h', col = cols, 
                 ylim = c(0, ymax)
                 , las = 1, xaxt='n', main =yvar)]
    abline(h=0, lty = 1)
        mtext(ll, outer = T, side = 1, line = -1)
    }
}
graphics.off()
## OK seems like I'm missing the m


with(itmp[lab==ll], plot(ordShow, spent_EV, xlab='individual', ylab='risk spent', bty = 'n', type = 'h', col = cols, ylim = ylim, las = 1, xaxt='n', main =ll, xlim = c(0,600)))

## risk averted & avertable risk not averted for each design
## avert_EV & spent_EV

## add error bars to inf spent & frac of information
## try fraction of information per person on y-axis
## show inf spent/averted on same plot

## equipoise perturbed plot
## histogram within strata groups
## look at risk by person/strata by treatment assignment for trials
## vaccinating a greater % of people increases risk averted/spent, but still has problem of withholding treatment from individuals

par(mfrow=c(3,1))
xmax <- 600
itmp[lab=='RCT',plot(ordShowArm.sp, inf, type = 'h', col = cols, xlim = c(0,xmax))]
itmp[lab=='RCT',plot(ordShowArm.sp, spent_EV, type = 'h', col = cols, xlim = c(0,xmax))]
itmp[lab=='RCT',plot(ordShowArm.sp, avert_EV, type = 'h', col = cols, xlim = c(0,xmax))]

