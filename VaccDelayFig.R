if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevenbellan', Sys.info()['user'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('sbellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('nid', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/EbolaVaccSim/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T)); gc()
require(optiRum); require(RColorBrewer); require(boot); require(data.table); require(ggplot2); require(grid); require(reshape2); require(parallel)
sapply(c('multiplot.R','extractFXN.R','ggplotTheme.R'), source)
sapply(c('simFuns.R','AnalysisFuns.R','ExpFit.R'), source)
avHazDo <- ''

####################################################################################################
## Delayed vaccination comparators & threshold-risk spending modifications for Supplementary Appendix
####################################################################################################
## Need to load both CVD and general simulations (careful not to double up on CVD which may be represeneted in both for some scenarios)
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
punq[,clusSize:=as.numeric(clusSize)]
plunq <- punq[threshold==.05, list(lab, power, trialStartDate, threshold, above, above_EV, caseSpent, totCase, totCase_EV,avHaz, tcal,date, clusSize)]
irsk <- merge(irsk, unique(punq[,.(tid, clusSize= as.numeric(clusSize))]), by='tid')## get clusSize in irsk for ordering purposes

numClus <- as.numeric(punq[,unique(numClus)])
hazT <- setClusHaz(makePop(makeParms(trialStartDate="2014-10-01",propInTrial=0.05,avHaz=avHazDo,indivRRSeed=7,HazTrajSeed=7, numClus=numClus)))$hazT

irsk[,ordX:=ordShowArm]
irsk[lab=='SWCT',ordX:=ordShow]

unique(parmsMat[,.(clusSize,avHaz, trialStartDate, .N), tid])
punq[threshold==.01,.(clusSize,avHaz, trialStartDate, lab,tid)]
names(irsk)

## Add a 6th of a cluster of space b/w each cluster for easier visualization
irsk[,ordX.sp:=ordX + floor((ordX-1)/clusSize)*(clusSize/6)]
avertableTab <- unique(irsk[lab %in% c('NT','VR'), .(arms = T, avertableRisk = inf[lab=='NT']-inf[lab=='VR'], ordX, ordX.sp, cVaccOrd, inf=inf[lab=='NT']), .(Oi,tid, clusSize)])
avertableTab2 <- merge(avertableTab[,.(Oi,tid,arms,avertableRisk,inf, clusSize, cVaccOrd)], unique(irsk[lab =='SWCT', .(Oi, ordX, ordX.sp), clusSize]), by =c('Oi', 'clusSize'))
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

tidDo <- 2
labDo <- c('RCT-gs-rp-cvd', 'RCT-gs-rp')
itmp <- irsk[clusSize==clusSizeDo & tid==tidDo & lab %in% labDo]
##tidDo <- itmp[,unique(tid)]
unique(itmp[,.(arm, armShown, type, lab, exmpl, tid)])
itmp <- itmp[tid==tidDo & ((arm==armShown & type=='cond' &  grepl('RCT',lab)) | (type=='condvd' & lab=='SWCT' & exmpl==T))]
itmp[, cols:=armShown]; itmp[cols=='cont',cols:='red']; itmp[cols=='vacc',cols:='dodger blue']
itmp <- itmp[lab!='RCT-rp']
unique(itmp[,.(arm, armShown, type, lab, exmpl, tid)])

atmp <- avertableTab[tid==tidDo & clusSize==clusSizeDo]


## for plotting
clusSizetmp <- punq[tid==tidDo, as.numeric(clusSize)][1]
Ntmp <- punq[tid==tidDo, as.numeric(clusSize)*as.numeric(numClus)][1]
mids <- seq(clusSizetmp/2, Ntmp-clusSizetmp/2, by = clusSizetmp)
mids <- mids + clusSizetmp/6*c(0:(length(mids)-1))
tcks <- rep(mids, each = 2) + c(-clusSizetmp/2,clusSizetmp/2)
ncol <- 2

wid <- 9
heig <- 6
res <- 300
opacity <- 100
mcex <- .7

labsToDo <- c('RCT-gs-rp','RCT-gs-rp-cvd')
lenL <- length(labsToDo)
labsDisplayRaw <-c( 'RCT (risk-prioritized rollout, three analyses)',
                   'RCT (risk-prioritized rollout, three analyses, 60-day delayed vacc. comparator')
labsDisplay <- paste0('(',LETTERS[1:lenL],') ', labsDisplayRaw)

 
numClus <- as.numeric(punq[,unique(numClus)])
lwdhbase <- 1.2
contcol <- 'red'
cvdcol <- 'purple'
####################################################################################################
## Effect of vacc delay by cluster size
for(jj in 1:2) {
    if(jj==2) {
        ylimavertable <- c(0, max(avertableTab[,avertableRisk]))
    }else{
        ylimavertable <- c(0, 1)
    }
    xlim <- c(0,max(mids) + 390)
    #pdf(file.path(figdir, paste0('effect of vacc delay', 'avertable'[jj==0], 'total'[jj==1], '.pdf')), w = wid, h = heig); lwdhbase <- lwdhbase/1.2
    png(file.path(figdir, paste0('Fig S4 VaccDel', 'avertable'[jj==0], 'total'[jj==1], '.png')), w = wid, h = heig, res = res, units = 'in')
    layout(matrix(c(1:5, 1:4, 6),5,2), heights = c(rep(1,4),3))
    par(mar = c(0,8,1,3), oma = c(1,1.5,2,.5), ps = 17)
    magnifier <- max(clusSizes)/clusSizes
    for(cc in 1:(length(clusSizes)-1)) {
        clusSizeDo <- clusSizes[cc]
        lwdh <- lwdhbase * (max(clusSizes)/clusSizeDo)^.06
        itmp <- irsk[clusSize==clusSizeDo & lab %in% labDo]
        tidDo <- punq[clusSize==clusSizeDo & avHaz==avHazDo & lab %in% labDo, unique(tid)]
        itmp <- itmp[tid==tidDo & ((arm==armShown & type=='cond' &  grepl('RCT',lab)) | (type=='condvd' & lab=='SWCT' & exmpl==T))]
        itmp[, cols:=armShown]; itmp[cols=='cont',cols:=contcol]; itmp[armShown=='cont' & grepl('cvd',lab) ,cols:=cvdcol]; itmp[cols=='vacc',cols:='dodger blue']
        itmp <- itmp[lab!='RCT-rp']
        atmp <- avertableTab[tid==tidDo & clusSize==clusSizeDo]
        if(jj==1) {
            irsk[lab=='NT' & clusSize==clusSizeDo & tid==tidDo, plot(magnifier[cc]*ordX.sp, inf, ylab ='', bty = 'n', type = 'h', lwd = lwdh, xaxt='n',
                                                                     lend=1, col = 'black', xlim = xlim, ylim = ylimavertable, main = , las = 1)]
            atmp[arms==T, points(magnifier[cc]*ordX.sp, avertableRisk, type = 'h', lwd = lwdh, lend=1, col = 'gray')]
        }else{
            atmp[arms==T, plot(magnifier[cc]*ordX.sp, avertableRisk, ylab ='', bty = 'n', type = 'h', lwd = lwdh, lend=1, col = 'gray', xlim = xlim, ylim = ylimavertable, axes=F,
                               main = '', las = 1)]
            axis(2, at = seq(0,.25, by = .05), labels = NA)
            axis(2, at = seq(0,.2, by = .1), lwd = 0,las=1)
        }
        itmp[lab=='RCT-gs-rp-cvd', points(magnifier[cc]*ordX.sp, avert_EV,  type = 'h', lwd = lwdh, lend=1, col = makeTransparent(cols,opacity))]
        itmp[lab=='RCT-gs-rp', points(magnifier[cc]*ordX.sp, avert_EV, type = 'h', lwd = lwdh, lend=1, col = makeTransparent(cols,opacity))]
        mtext(clusSizeDo, side = 2, 6, las = 1, cex=mcex)
        xbump <- 160
        xbump2 <- 300
        ## power
        points(max(tcks)+xbump+25, plunq[avHaz==avHazDo & clusSize==clusSizeDo & lab=='RCT-gs-rp', power]*ylimavertable[2], type = 'h', lwd = lwdh*7, lend=1, col = contcol)
        points(max(tcks)+xbump+65, plunq[avHaz==avHazDo & clusSize==clusSizeDo & lab=='RCT-gs-rp-cvd', power]*ylimavertable[2], type = 'h', lwd = lwdh*7, lend=1, col = cvdcol)
        ## duration
        points(max(tcks)+xbump2+25, plunq[avHaz==avHazDo & clusSize==clusSizeDo & lab=='RCT-gs-rp', tcal]/168*ylimavertable[2], type = 'h', lwd = lwdh*7, lend=1, col = contcol)
        points(max(tcks)+xbump2+65, plunq[avHaz==avHazDo & clusSize==clusSizeDo & lab=='RCT-gs-rp-cvd', tcal]/168*ylimavertable[2], type = 'h', lwd = lwdh*7, lend=1, col = cvdcol)
        axis(4, at = ylimavertable[2]*seq(0,1,by=.5), labels=seq(0,1,by=.5)*24, line = .3, las = 1)
        axis(2, at = ylimavertable[2]*seq(0,1,by=.5), labels=seq(0,1,by=.5), line = -50.5, las = 1)
        mtext(paste0('(',LETTERS[cc],')'), side = 3, line = 0, adj = .01, cex=mcex) 
        for(ii in 1:numClus) segments(tcks[2*ii-1],0, tcks[2*ii],0, lwd = .7)
        par(xpd=T)
        if(cc==1) legend(800, ifelse(jj==2, .35, 1.3), c('avertable', 'averted (vaccine arm)','averted (control)', 'averted (60-day delayed vacc. comparator)'), 
                         col = c('gray','blue','red','purple'), bty = 'n', pch = 15, cex =  1)
        par(xpd=F)
    }
    axis(1, at = max(tcks) + xbump + 37.5, lab = c('power'), las = 2, lwd = 0)
    axis(1, at = max(tcks) + xbump2 + 37.5, lab = c('duration\n(weeks)'), las = 2, lwd = 0)#, col.axis = 'dark green')
    axis(1, at = mids, lab = 1:numClus, lwd = 0)
    mtext(paste0('individuals (grouped by cluster and arm)'), 1, outer=F, line = 3, cex=mcex)
    mtext('individual risk', 2,outer=T, line = -4, cex=mcex, adj = .75)
    par(xpd=T)
    mtext('cluster\nsize', 3, outer = T, adj = 0, line = -1, cex = mcex)
    plunq[,col:='black']
    plunq[lab==labsToDo[1],col:='red']
    plunq[lab==labsToDo[2],col:='purple']
    ## power vs clusSize
    par(mar = c(4,5,7,1))
    plunq[avHaz== avHazDo & lab %in% labsToDo, plot(clusSize, power, col = col, type = 'n', xlab = '', ylab = 'power', bty = 'n',las=1, axes=F, ylim = c(.5,.9))]
    plunq[avHaz== avHazDo & lab %in% labsToDo, points(clusSize, power, col = col, type = 'b', pch = 16), lab]
    axis(1, clusSizes)
    axis(2, c(.5,.7,.9), las = 1)
    legend('bottomright', c('control arm', '60-day delayed vacc. comparator arm'), col = c('red','purple'), bty = 'n', pch = 16, cex =  1)
    mtext('(E)', side = 3, line = 0, adj = .01, cex=mcex) 
    ## speed vs clusSize
    plunq[avHaz== avHazDo & lab %in% labsToDo, plot(clusSize, tcal/7, col = col, type = 'n', xlab = '', ylab = 'duration (weeks)', bty = 'n',las=1, axes=F, ylim = c(24,0))]
    plunq[avHaz== avHazDo & lab %in% labsToDo, points(clusSize, tcal/7, col = col, type = 'b', pch = 16), lab]
    axis(1, clusSizes)
    axis(2, at = seq(0,24,4), las = 1)
    mtext('(F)', side = 3, line = 0, adj = .01, cex=mcex) 
    par(xpd=T)
    arrows(205,0,205,24, len = .1, code = 1)
    mtext('speed', 4, cex=mcex)
    mtext('individuals per cluster', 1, outer=T, line = -.5, cex = mcex)
    graphics.off()
 }

plunq[avHaz==avHazDo &lab %in% labsToDo, .(clusSize, lab, power, tcal)]

## standalone legend
pdf(file.path(figdir, paste0('vaccDel leg.pdf')), width = 3, height = 2)
plot(0,0,type='n', bty ='n', axes =F, xlab='',ylab='')
par(xpd=T)
legend('center', col = c('gray', contcol, cvdcol, 'dodger blue'), 
       leg = c('avertable', 'averted (control)', 'averted (60-day delayed vaccine)', 'averted (vaccine)'), pch = 15, cex = .8, bty = 'n')
dev.off()

## standalone legend
pdf(file.path(figdir, paste0('vaccDel leg total.pdf')), width = 3, height = 2)
plot(0,0,type='n', bty ='n', axes =F, xlab='',ylab='')
par(xpd=T)
legend('topright', col = c('black','gray', contcol, cvdcol, 'dodger blue'), 
       leg = c('total', 'avertable', 'averted (control)', 'averted (60-day delayed vaccine)', 'averted (vaccine)'), pch = 15, cex = .8, bty = 'n')
dev.off()
 



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
        lwdh <- lwdhbase * max(clusSizes)/clusSizeDo
        itmp <- irsk[clusSize==clusSizeDo]
        tidDo <- itmp[,unique(tid)]
        itmp <- itmp[tid==tidDo & ((arm==armShown & type=='cond' &  grepl('RCT',lab)) | (type=='condvd' & lab=='SWCT' & exmpl==T))]
        itmp[, cols:=armShown]; itmp[cols=='cont',cols:='red']; itmp[cols=='vacc',cols:='dodger blue']
        itmp <- itmp[lab!='RCT-rp']
        atmp <- avertableTab[tid==tidDo & clusSize==clusSizeDo]
        if(jj==1) {
            irsk[lab=='NT' & clusSize==clusSizeDo & tid==tidDo, plot(magnifier[cc]*ordX.sp, inf, ylab ='', bty = 'n', type = 'h', lwd = lwdh, lend=1, col = 'black', xlim = xlim, ylim = ylimavertable, xaxt='n', 
                                                                     main = paste0(clusSizeDo, '/cluster'), las = 1)]
            atmp[arms==T, points(magnifier[cc]*ordX.sp, avertableRisk, type = 'h', lwd = lwdh, lend=1, col = 'gray')]
        }else{
            atmp[arms==T, plot(magnifier[cc]*ordX.sp, avertableRisk, ylab ='', bty = 'n', type = 'h', lwd = lwdh, lend=1, col = 'gray', xlim = xlim, ylim = ylimavertable, xaxt='n', 
                               main = paste0(clusSizeDo, '/cluster'), las = 1)]
        }
#        itmp[lab=='RCT-gs-rp-cvd', points(magnifier[cc]*ordX.sp, avert_EV,  type = 'h', lwd = lwdh, lend=1, col = makeTransparent(cols,opacity))]
        itmp[lab=='RCT-gs-rp', points(magnifier[cc]*ordX.sp, avert_EV, type = 'h', lwd = lwdh, lend=1, col = makeTransparent(cols,opacity*2))]
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
