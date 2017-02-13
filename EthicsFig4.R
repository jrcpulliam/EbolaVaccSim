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

thing <- 'Equip-Fig5-v2'
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
plunq <- punq[threshold==.05 & lab %in% labstoshow, list(lab, power, trialStartDate, threshold, above, above_EV, caseSpent_EV, totCase, totCase_EV,avHaz, tcal,date, cols, pchs)]

## ## ethicline contours
## ias <- pretty(plunq$caseSpent, n = 50)
## pows <- seq(0,1,l=50)
## contt <- data.table(expand.grid(ia=ias, pow=pows))
## ii=1
## infPpow <- c(100, 200, 300)[ii]
## contt[, merit:= -ia/infPpow + pow] ## every case spent is worth 10% power
## ##    contt[pow<.3, merit:= -ia/infPpow+.3] ## don't penalize less power under 30% power since it's so insignificant anyways
## tmp <-  plunq[lab!='RCT-rp' & avHaz=='' & date=='Dec-14']
## ## totcase vs power
## p <- ggplot(tmp) + 
##     geom_tile(data = contt, aes(x=ia, y=pow, z=merit, fill=merit)) + scale_fill_gradient(low = "light green", high = "dark green") +
##         stat_contour(data = contt, aes(x=ia, y=pow, z=merit), col = gray(.9, alpha = .9), bins = 7) + 
##             geom_point(data = tmp[date!='Apr-15'], aes(caseSpent, power, colour = lab, shape = lab), size =2.5) +
##                 ylim(0,1) + xlab(paste0('cases spent relative to vaccine rollout')) +
##                     scale_color_manual(values=cols) + scale_shape_manual(values=pchs) + 
##                         ylab('power \n(given efficacy>0)') + theme(axis.title.y = element_text(angle=0))
## ggsave(file.path(figdir, paste0('power vs case.pdf')), plot=p, w=wid, h=heig, units='in', dpi = 200)


## Fig 4B with just power vs avertable cases that aren't averted
tmpM <-  plunq[lab!='RCT-rp' & avHaz=='' & date %in% c('Dec-14','Feb-15')]
## tmp[date=='Dec-14', cols:=makeTransparent(cols, 100)]
##tmp[date=='Feb-15',pch:=15]
## tmpM[pchs==21 & date=='Dec-14', pchs:=16]
## tmpM[pchs==21 & date=='Dec-14', pchs:=15]
## tmpM[pchs==21 & date=='Dec-14', pchs:=17]

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

pdf(file.path(figdir, paste0('tst.pdf')), width =5.5, h=3)
plot(1:30, pch = 1:30)
graphics.off()

load(file='data/WAevddat.Rdata')
SLdistr[col=

SLdistr[, plot(date, cases, type = 'n', las = 1)]
SLdistr[,lines(date, cases, col = Location)]

pdf(file.path(figdir, paste0('Fig4B.pdf')), width =6.5, h=3)
par(mfrow=c(1,3), ps='14')
par(mar = c(4 ,6,1,.5))#par(mar=c(5,5,.5,.5))
hazT$lwd <- 1
hazT[,col:=makeTransparent(rainbow(length(unique(cluster)))[cluster], alpha = 180)]
hazT[,plot(Date,clusHaz, lwd=0, bty = 'n', las = 1, ylab = '', type = 'n', xlab='')]
hazT[, lines(Date, clusHaz, lwd = lwd, col = col), cluster]
par(xpd=T)
rect(head(tmpM$trialStartDate,1), -.00095, head(tmpM$trialStartDate,1)+24*7, -.0008, col = 'purple', border=NA)
rect(tail(tmpM$trialStartDate,1), -.0010, tail(tmpM$trialStartDate,1)+24*7, -.0012, col = 'orange', border=NA)
#mtext('daily hazard', 2, 4, cex = .7)
title(ylab='daily hazard', line = 4)
title(main = 'A')
cex <- 2
cexleg <- .8
lwd <- 2
par(mar = c(4 ,4,1,.5))                
tmpM[,pchs2:=as.numeric(as.character(factor(lab, labels=c(2,6:5,3:4,16))))]
tmpM[,cols2:=makeTransparent(ifelse(date==date[1], 'purple','orange'), alpha = 160)]
xmax <- 45
## power
plot(0, 0, type = 'n', xlim = c(0,xmax), ylim = c(0,1), las = 1, bty = 'n', ylab = 'power', xlab = '')
tmpM[, points(caseSpent_EV, power, xlim = c(0,xmax), ylim = c(0,1), pch = pchs2, las = 1, col = cols2, cex = cex, lwd=lwd), date]
tmpM[date==date[1], legend('bottomright', leg = lab, pch = pchs2, cex = cexleg, bty = 'n')]
tmpM[lab==lab[1], legend('bottomleft', leg = date, pch = 16, col=cols2, cex = cexleg, bty = 'n', title = 'start date')]
title(main = 'B')
## speed
plot(0, 0, type = 'n', xlim = c(0,xmax), ylim = c(24,0), las = 1, bty = 'n', ylab = 'trial duration (weeks)', xlab = '', yaxt = 'n')
axis(2, at = 0:6*4, las = 1)
tmpM[, points(caseSpent_EV, tcal/7, pch = pchs2, col = cols2, cex = cex, lwd=lwd), date]
mtext('cases spent relative to optimized rollout', 1, -1, outer = T, at = .65, cex = .7)
title(main = 'C')
graphics.off()
  

 

