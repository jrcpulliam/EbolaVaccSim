if(grepl('Stevens-MacBook', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevenbellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('nid', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/EbolaVaccSim/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T)); gc()
require(RColorBrewer); require(boot); require(data.table); require(ggplot2); require(grid); require(reshape2); require(parallel); require(optiRum); 
sapply(c('multiplot.R','extractFXN.R','ggplotTheme.R'), source)
wid <- 6.5
heig <- 4
res <- 300

## think about cRCT and rpSWCT on this plot. otherwise kinda finished??
sapply(c('simFuns.R','AnalysisFuns.R','ExpFit.R'), source)
hazT <- setClusHaz(makePop(makeParms(trialStartDate="2014-10-01",propInTrial=0.05,avHaz="",indivRRSeed=7,HazTrajSeed=7)))$hazT

thing <- 'Equip-Fig5-10clus-3cat'
figdir <- file.path('Figures', thing)
dir.create(figdir)
fls <- list.files('BigResults', pattern = paste0(thing,'-'), full.names=T)

punq <- data.table()
for(ii in 1:length(fls)) {
    load(fls[ii])
    tmp <- merge(resList$punq,resList$SpopWorst,by = c('pid', 'propInTrial','lab'))
    tmp <- merge(tmp, resList$irsk[type=='marg', list(caseSpent_EV= 6000*mean(spent_EV), caseSpent= 6000*mean(spent)), pid], by = 'pid')
    speedpow <- with(resList, merge(parms[, list(pid, nbatch)], finTrials[,list(tcal, vaccEff,vaccGood, cvr,nbatch)], by = 'nbatch'))
    tmp <- merge(tmp, speedpow[, list(tcal = mean(tcal)), pid], by = 'pid')
    punq <- rbind(punq, tmp)
}
punq[,trialStartDate:=as.Date(trialStartDate)]
punq[,date:=format.Date(trialStartDate, '%b-%y')]

cols <- c("NT"='red', "SWCT"='orange', "VR"='dark green', "RCT-gs"='dodger blue', "RCT"='light blue', "RCT-rp"='purple', "RCT-gs-rp" = 'purple', 'RCT-gs-rp-cvd'='purple', 'RCT-gs-rp-maxRR'='purple')
punq$cols <- cols[as.character(punq[,lab])]
punq$pchs <- 21
punq[lab=='RCT-gs-rp-cvd',pchs:=22]
punq[lab=='RCT-gs-rp-maxRR',pchs:=24]
## Fig 4
labstoshow <- c('RCT' ,'RCT-gs', 'RCT-rp', 'RCT-gs-rp', 'SWCT') 
labstoshow <- punq[,unique(lab)]
plunq <- punq[threshold==.05 & lab %in% labstoshow, list(lab, power, trialStartDate, threshold, above, above_EV, caseSpent_EV, totCase, totCase_EV,avHaz, tcal,date, cols, pchs, clusSize)]

names(punq)
plunq

## Fig 4B with just power vs avertable cases that aren't averted
tmpM <-  plunq[lab!='RCT-rp' & avHaz==''  & clusSize==200]

cex <- 2
#pdf(file.path(figdir, paste0('Fig4B.pdf')), width =5.5, h=3)
par(mfrow=c(1,2), mar = c(4 ,6,.5,.5))
par(mfrow=c(2,2), mar = c(3 ,6,.5,.5), oma = c(1,0,0,0))
for(dd in tmpM[,unique(date)]) {
    tmp <- tmpM[date==dd]
    xmax <- 60
    ## power
    with(tmp, plot(caseSpent_EV, power, xlim = c(0,xmax), ylim = c(0,1), pch = pchs, las = 1, col = cols, cex = cex, bty = 'n', xlab = '', main = date[1]))
    mtext('avertable cases not averted', 1, -1, outer = T)
    ## speed
    with(tmp, plot(caseSpent_EV, tcal/30, xlim = c(0,xmax), ylim = c(180/30,0), pch = pchs, las = 1, col = cols, cex=cex, bty='n', xlab='', ylab = 'trial duration\n(months)'))
}
    tmp[, legend('bottomright', leg = lab, col = cols, pch = pchs, cex = 1, bty = 'n')]
#graphics.off()

load(file='data/WAevddat.Rdata')

tmpM[,pchs2:=as.numeric(as.character(factor(lab, labels=c(2,6:5,3:4,16))))]
tmpM[lab=='RCT', pchs2:=13]
tmpM[,cols2:=makeTransparent(ifelse(date=='Oct-14', 'purple','orange'), alpha = 160)]

mtextadj <- -.3
mline <- .7
pdf(file.path(figdir, paste0('Fig4.pdf')), width =5.5, 5)
par(mfrow=c(1,3), ps='14')
layout(matrix(c(1,2,4,3),2,2))
par(mar = c(3.5 ,5,1.5,0), oma = c(0,0,.5,2))
hazT$lwd <- 1
hazT[,col:=makeTransparent(rainbow(length(unique(cluster)))[cluster], alpha = 180)]
hazT[,plot(Date,clusHaz, lwd=0, bty = 'n', las = 1, ylab = '', type = 'n', xlab='')]
hazT[, lines(Date, clusHaz, lwd = lwd, col = col), cluster]
par(xpd=T)
rect(tail(tmpM$trialStartDate,1), -.0013 , tail(tmpM$trialStartDate,1)+24*7, -.0015, col = 'purple', border=NA)
rect(head(tmpM$trialStartDate,1), -.00155, head(tmpM$trialStartDate,1)+24*7, -.00175, col = 'orange', border=NA)
#mtext('daily hazard', 2, 4, cex = .7)
title(ylab='daily hazard', line = 4)
mtext('(A)', side = 3, line = mline, adj = mtextadj)
cex <- 2
mcex <- .7
cexleg <- .8
lwd <- 2
#par(mar = c(4 ,4,1,1))
xmax <- 210
## power 
plot(0, 0, type = 'n', xlim = c(0,xmax), ylim = c(0,1), las = 1, bty = 'n', ylab = '', xlab = '', xaxt='n')
title(ylab='power', line = 4)
axis(1, seq(0,200, by = 50))
tmpM[, points(caseSpent_EV, power, xlim = c(0,xmax), ylim = c(0,1), pch = pchs2, las = 1, col = cols2, cex = cex, lwd=lwd), date]
#tmpM[date==date[1], legend('bottomright', leg = lab, pch = pchs2, cex = cexleg, bty = 'n')]
#tmpM[lab==lab[1], legend('bottomleft', leg = date, pch = 16, col=cols2, cex = cexleg, bty = 'n', title = 'start date')]
mtext('(C)', side = 3, line = mline, adj = mtextadj)
## speed
plot(0, 0, type = 'n', xlim = c(0,xmax), ylim = c(24,0), las = 1, bty = 'n', ylab = 'trial duration (weeks)', xlab = '', axes=F)
axis(1, seq(0,200, by = 50))
axis(2, at = 0:6*4, las = 1)
tmpM[, points(caseSpent_EV, tcal/7, pch = pchs2, col = cols2, cex = cex, lwd=lwd), date]
mtext('cases spent relative to optimized rollout', 1, -1.5, outer = T, at = .5, cex = mcex)
par(xpd=NA)
arrows(xmax*1.05, 0, xmax*1.05, 24, code = 1, len = .1)
mtext('trial speed', 4, .5, cex = .85)
mtext('(D)', side = 3, line = mline, adj = mtextadj)
plot(0, 0, type = 'n', xlab='',ylab='', bty='n', axes=F)
mtext('(B)', side = 3, line = mline, adj = mtextadj)
graphics.off()
     
save(tmpM, file = 'BigResults/Fig4data.Rdata')
 

