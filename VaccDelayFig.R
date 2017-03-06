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
hazT <- setClusHaz(makePop(makeParms(trialStartDate="2014-10-01",propInTrial=0.05,avHaz="xTime",indivRRSeed=7,HazTrajSeed=7)))$hazT

thing <- 'Equip-Fig-SX-vaccDel'
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

irsk[,ordX:=ordShowArm]
irsk[lab=='SWCT',ordX:=ordShow]

unique(parmsMat[,.(clusSize,avHaz, trialStartDate, .N), tid])
punq[threshold==.01,.(clusSize,avHaz, trialStartDate, lab,tid)]
names(irsk)

## Add a 6th of a cluster of space b/w each cluster for easier visualization
irsk[,ordX.sp:=ordX + floor((ordX-1)/clusSize)*(clusSize/6)]
avertableTab <- unique(irsk[lab %in% c('NT','VR'), .(arms = T, avertableRisk = inf[lab=='NT']-inf[lab=='VR'], ordX, ordX.sp, cVaccOrd, inf=inf[lab=='NT']), .(Oi,tid, clusSize)])
avertableTab2 <- merge(avertableTab[,.(Oi,tid,arms,avertableRisk,inf, cVaccOrd)], unique(irsk[lab =='SWCT', .(Oi, ordX, ordX.sp)]), by ='Oi')
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


clusSizeDo <- 200                 
if(!clusSizeDo %in% clusSizes) clusSizeDo <- clusSizes[1]
cc <- which(clusSizes==clusSizeDo)
## unique(irsk[,.(arm, armShown, type, lab, exmpl, clusSize,tid)])

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

numClus <- as.numeric(punq[,unique(numClus)])
####################################################################################################
## Effect of vacc delay by cluster size
for(jj in 1:2) {
    if(jj==2) {
        ylimavertable <- c(0, max(avertableTab[,avertableRisk]))
    }else{
        ylimavertable <- c(0, 1)
    }
    xlim <- c(0,Ntmp)
#    pdf(file.path(figdir, paste0('effect of vacc delay', 'avertable'[jj==0], 'total'[jj==1], '.pdf')), w = wid, h = heig)
    png(file.path(figdir, paste0('effect of vacc delay', 'avertable'[jj==0], 'total'[jj==1], '.png')), w = wid*100, h = heig*100)
    par(mfrow=c(4,1), mar = c(0,5,1,0), oma = c(5,2,0,0), ps = 17)
    magnifier <- max(clusSizes)/clusSizes
    for(cc in 1:(length(clusSizes)-1)) {
        clusSizeDo <- clusSizes[cc]
        itmp <- irsk[clusSize==clusSizeDo]
        tidDo <- itmp[,unique(tid)]
        itmp <- itmp[tid==tidDo & ((arm==armShown & type=='cond' &  grepl('RCT',lab)) | (type=='condvd' & lab=='SWCT' & exmpl==T))]
        itmp[, cols:=armShown]; itmp[cols=='cont',cols:='red']; itmp[cols=='vacc',cols:='dodger blue']
        itmp <- itmp[lab!='RCT-rp']
        atmp <- avertableTab[tid==tidDo & clusSize==clusSizeDo]
        if(jj==1) {
            irsk[lab=='NT' & clusSize==clusSizeDo & tid==tidDo, plot(magnifier[cc]*ordX.sp, inf, ylab ='', bty = 'n', type = 'h', col = 'black', xlim = xlim, ylim = ylimavertable, xaxt='n', 
                                                                     main = paste0(clusSizeDo, '/cluster'), las = 1)]
            atmp[arms==T, points(magnifier[cc]*ordX.sp, avertableRisk, type = 'h', col = 'gray')]
        }else{
            atmp[arms==T, plot(magnifier[cc]*ordX.sp, avertableRisk, ylab ='', bty = 'n', type = 'h', col = 'gray', xlim = xlim, ylim = ylimavertable, xaxt='n', 
                               main = paste0(clusSizeDo, '/cluster'), las = 1)]
        }
        itmp[lab=='RCT-gs-rp-cvd', points(magnifier[cc]*ordX.sp, avert_EV,  type = 'h', col = makeTransparent(cols,opacity))]
        itmp[lab=='RCT-gs-rp', points(magnifier[cc]*ordX.sp, avert_EV, type = 'h', col = makeTransparent(cols,opacity*2))]
    }
    if(jj==2) {
legend('topright', col = c('gray', makeTransparent("red", c(opacity*2, opacity)), 'dodger blue'),
           leg = c('avertable', 'averted (control)', 'averted (60-day delayed vaccine)', 'averted (vaccine)'), pch = 15, cex = 1.3, bty = 'n')
}else{
    legend('topright', col = c('black','gray', makeTransparent("red", c(opacity*2, opacity)), 'dodger blue'),
           leg = c('total', 'avertable', 'averted (control)', 'averted (60-day delayed vaccine)', 'averted (vaccine)'), pch = 15, cex = 1.3, bty = 'n')
}
    axis(1, at = mids, lab = 1:numClus, lwd = 0)
    mtext(paste0('individuals by cluster'), 1, outer=T, line = 3)
    mtext('risk', 2, ,outer=T, line = -1)
    graphics.off()
}


## For PPT

numClus <- as.numeric(punq[,unique(numClus)])
####################################################################################################
## Effect of vacc delay by cluster size
for(jj in 1:2) {
    if(jj==2) {
        ylimavertable <- c(0, max(avertableTab[,avertableRisk]))
    }else{
        ylimavertable <- c(0, 1)
    }
    xlim <- c(0,Ntmp)*1.2
    pdf(file.path(figdir, paste0('effect of vacc delay', 'avertable'[jj==0], 'total'[jj==1], '.pdf')), w = wid, h = heig)
#    png(file.path(figdir, paste0('effect of vacc delay', 'avertable'[jj==0], 'total'[jj==1], '.png')), w = wid*100, h = heig*100)
    par(mfrow=c(1,1), mar = c(0,5,1,1), oma = c(5,2,0,0), ps = 17)
    magnifier <- max(clusSizes)/clusSizes
    for(cc in 1) {
        clusSizeDo <- clusSizes[cc]
        itmp <- irsk[clusSize==clusSizeDo]
        tidDo <- itmp[,unique(tid)]
        itmp <- itmp[tid==tidDo & ((arm==armShown & type=='cond' &  grepl('RCT',lab)) | (type=='condvd' & lab=='SWCT' & exmpl==T))]
        itmp[, cols:=armShown]; itmp[cols=='cont',cols:='red']; itmp[cols=='vacc',cols:='dodger blue']
        itmp <- itmp[lab!='RCT-rp']
        atmp <- avertableTab[tid==tidDo & clusSize==clusSizeDo]
        if(jj==1) {
            irsk[lab=='NT' & clusSize==clusSizeDo & tid==tidDo, plot(magnifier[cc]*ordX.sp, inf, ylab ='', bty = 'n', type = 'h', col = 'black', xlim = xlim, ylim = ylimavertable, xaxt='n', 
                                                                     main = paste0(clusSizeDo, '/cluster'), las = 1)]
            atmp[arms==T, points(magnifier[cc]*ordX.sp, avertableRisk, type = 'h', col = 'gray')]
        }else{
            atmp[arms==T, plot(magnifier[cc]*ordX.sp, avertableRisk, ylab ='', bty = 'n', type = 'h', col = 'gray', xlim = xlim, ylim = ylimavertable, xaxt='n', 
                               main = paste0(clusSizeDo, '/cluster'), las = 1)]
        }
#        itmp[lab=='RCT-gs-rp-cvd', points(magnifier[cc]*ordX.sp, avert_EV,  type = 'h', col = makeTransparent(cols,opacity))]
        itmp[lab=='RCT-gs-rp', points(magnifier[cc]*ordX.sp, avert_EV, type = 'h', col = makeTransparent(cols,opacity*2))]
    }
    if(jj==2) {
legend('topright', col = c('gray', makeTransparent("red", c(opacity*2, opacity)), 'dodger blue'),
           leg = c('avertable', 'averted (control)', 'averted (60-day delayed vaccine)', 'averted (vaccine)'), pch = 15, cex = 1.3, bty = 'n')
}else{
    legend('topright', col = c('black','gray', makeTransparent("red", c(opacity*2, opacity)), 'dodger blue'),
           leg = c('total', 'avertable', 'averted (control)', 'averted (60-day delayed vaccine)', 'averted (vaccine)'), pch = 15, cex = 1.3, bty = 'n')
}
    axis(1, at = mids, lab = 1:numClus, lwd = 0)
    mtext(paste0('individuals by cluster'), 1, outer=T, line = 3)
    mtext('risk', 2, ,outer=T, line = -1)
    graphics.off()
}


## load('data/vaccProp1.Rdata')
## vaccProp1[,mean(vaccEff)]
## vaccProp1[vaccEff!=0,mean(vaccEff)]
## vaccProp1[, mean(vaccEff==0)]
## hist(vaccProp1$vaccEff)
