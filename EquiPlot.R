
if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('nid', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/EbolaVaccSim/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T)); gc()
require(optiRum); require(RColorBrewer); require(boot); require(data.table); require(ggplot2); require(grid); require(reshape2); require(parallel)
sapply(c('multiplot.R','extractFXN.R','ggplotTheme.R'), source)
wid <- 6.5
heig <- 4
res <- 300

thing <- 'Equip-RRcat'
figdir <- file.path('Figures', thing)
dir.create(figdir)
fls <- list.files('BigResults', pattern = paste0(thing,'-'), full.names=T)
cols <- c("NT"='red', "SWCT"='orange', "VR"='dark green', "RCT-gs"='dodger blue', "RCT"='light blue', "RCT-rp"='purple', "RCT-gs-rp" = 'purple')

punq <- data.table()
for(ii in 1:length(fls)) {
    load(fls[ii])
    tmp <- merge(resList$punq,resList$SpopWorst,by = c('pid', 'propInTrial','lab'))
    tmp <- merge(tmp, resList$irsk[type=='marg', list(caseSpent = 6000*mean(spent_EV)), pid], by = 'pid')
    speedpow <- with(resList, merge(parms[, list(pid, nbatch)], finTrials[,list(tcal, vaccEff,vaccGood, cvr,nbatch)], by = 'nbatch'))
    tmp <- merge(tmp, speedpow[, list(tcal = mean(tcal)), pid], by = 'pid')
    punq <- rbind(punq, tmp)
}
punq[,trialStartDate:=as.Date(trialStartDate)]
punq[,date:=format.Date(trialStartDate, '%b-%y')]
plunq <- punq[threshold==.05, list(lab, power, trialStartDate, threshold, above, above_EV, caseSpent, totCase, totCase_EV,avHaz, tcal,date)]


load(file=paste0('BigResults/Equip-indivL/hazT7.Rdata')) ## make changeable later
setkey(hazT, day, cluster)
rect <- data.table(xmin=as.Date(punq[,unique(trialStartDate)]), ymin=-Inf, ymax=.001*(5:2))
rect[,xmax :=xmin+24*7]
hazT$date <- format.Date(hazT$Date,"%b-%Y")
ip0 <- ggplot(hazT, aes(x=Date, y=clusHaz, col=as.factor(cluster))) + geom_line() + eb +
    scale_x_date(limits=as.Date(c('2014-08-01','2015-10-01')), date_labels = "%b-%y") + labs(ylab='relative hazard')
ip <- ip0 + eb + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey20", alpha=0.2, inherit.aes = FALSE)
ggsave(file.path(figdir, paste0('hazT.jpeg')), ip, w=4, h = 2)

## ethicline contours
ias <- pretty(plunq$caseSpent, n = 50)
pows <- seq(0,1,l=50)
contt <- data.table(expand.grid(ia=ias, pow=pows))
ii=1
infPpow <- c(100, 200, 300)[ii]
contt[, merit:= -ia/infPpow + pow] ## every case spent is worth 10% power
##    contt[pow<.3, merit:= -ia/infPpow+.3] ## don't penalize less power under 30% power since it's so insignificant anyways
tmp <-  plunq[lab!='RCT-rp' & avHaz=='']
## totcase vs power
p <- ggplot(tmp[date!='Apr-15']) +
    geom_tile(data = contt, aes(x=ia, y=pow, z=merit, fill=merit)) + scale_fill_gradient(low = "light green", high = "dark green") +
        stat_contour(data = contt, aes(x=ia, y=pow, z=merit), col = gray(.9, alpha = .9), bins = 7) + 
            geom_point(data = tmp[date!='Apr-15'], aes(caseSpent, power, colour = lab), size =2.5) +
                ylim(0,1) + xlab(paste0('cases spent relative to vaccine rollout')) +
                facet_grid(trialStartDate~.) + scale_color_manual(values=cols) +
                    ylab('power \n(given efficacy>0)') + theme(axis.title.y = element_text(angle=0))
ggsave(file.path(figdir, paste0('power vs case.pdf')), plot=p, w=wid, h=heig, units='in', dpi = 200)



p1 <- ggplot(tmp, aes(trialStartDate, power, col = lab)) +
     scale_color_manual(values=cols, guide=F) + #scale_x_date(date_labels = '%b-%y') +
        geom_line() + ylim(0,1) + xlab('')
p2 <- ggplot(tmp, aes(trialStartDate, above, col = lab)) +
     scale_color_manual(values=cols, guide=F) + #scale_x_date(date_labels = '%b-%y') +
        geom_line()  + xlab('') + ylab('max > .05 risk')
p3 <- ggplot(tmp, aes(trialStartDate, tcal, col = lab)) +
     scale_color_manual(values=cols, guide=F) + #scale_x_date(date_labels = '%b-%y') +
        geom_line()  + xlab('') + ylab('days') + ylim(0,168)
p4 <- ggplot(tmp, aes(trialStartDate, caseSpent, col = lab)) +
     scale_color_manual(values=cols, guide=F) + #scale_x_date(date_labels = '%b-%y') +
        geom_line()  + xlab('') + ylab('excess cases')
pdf(file.path(figdir, paste0('speed power above by time.pdf')), w = wid, h=heig)#, units='in', res = 300)
multiplot(p1,p2,p3, p4, layout = matrix(1:4, ncol=2))
graphics.off()


p <- ggplot(tmp, aes(above, power, col = lab)) +
    facet_grid(date~.) + scale_color_manual(values=cols) +
        geom_point() + ylim(0,1) + xlab(paste0('maximum # subjects spending > threshold infection risk'))
ggsave(file.path(figdir, paste0('pow trhes.pdf')), plot=p, w=wid, h=heig, units='in')

p <- ggplot(plunq[lab!='RCT-rp'], aes(trialStartDate, power, col = lab, linetype=avHaz)) +
        geom_line() + ylim(0,1) + xlab('')
ggsave(file.path(figdir, paste0('pow trhes.pdf')), plot=p, w=wid, h=heig, units='in')

p <- ggplot(punq[lab!='RCT-rp'], aes(trialStartDate, above, col = lab, linetype=avHaz)) +
        geom_line() + xlab('') + facet_grid(threshold~.)
ggsave(file.path(figdir, paste0('pow trhes.pdf')), plot=p, w=wid, h=heig, units='in')

load(fls[grepl(23, fls)])
attach(resList)

####################################################################################################
####################################################################################################
lbrks <- c(.0001,.0005,.001,.005,.01,.05,.1,.5,1)

## cols <- c("NT"='red', "SWCT"='orange', "VR"='dark green', "RCT-gs"='dodger blue', "RCT"='blue', "RCT-rp"='purple', "RCT-gs-rp" = 'magenta')

jpeg(file.path(figdir, paste0('Hpop.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(Hpop[mids>=.005], aes(mids, counts, col=lab)) + geom_line() + xlim(0, .2) + ylim(0,2000)
graphics.off()

jpeg(file.path(figdir, paste0('Hpop_EV.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(Hpop[mids>=.005], aes(mids, counts_EV, col=lab)) + geom_line() + xlim(0, .2) + ylim(0,2000)
graphics.off()

jpeg(file.path(figdir, paste0('Hpop by EV.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(Hpop[mids>=.005]) + geom_line(aes(mids, counts, col=lab)) +
    geom_line(aes(mids, counts_EV, col=lab), linetype = 2) + xlim(.02-.005, .2) + ylim(0,500)
graphics.off()

## density lines risk spent
tmp <- irsk[type=='marg' & !lab%in%c('NT','VR')]
adj <- 7
p <- ggplot() + 
    geom_line(data=tmp, aes(x=spent_EV, col=lab, linetype = gs), stat='density', adjust=adj) +
              scale_color_manual(values=cols)  + xlab('per capita infection risk spent')  + ggtitle('average risk spent per individual') #+ xlim(.005, .3) + ylim(0,32)
ggsave(file.path(figdir, paste0('irsk spentEV.jpeg')), p, w = wid, h = heig, units = 'in')
