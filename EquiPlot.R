if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('nid', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/EbolaVaccSim/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T)); gc()
require(optiRum); require(RColorBrewer); require(boot); require(data.table); require(ggplot2); require(grid); require(reshape2); require(parallel)

## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('multiplot.R','extractFXN.R','ggplotTheme.R'), source)
wid <- 6.5
heig <- 4
res <- 300
## Make Figure Folder
figdir <- file.path('Figures', thing)
dir.create(figdir)

## thing <- 'Equip-RRcat'
thing <- 'Equip-irskHybrid'
load(file=file.path('BigResults',paste0(thing, '.Rdata')))
names(resList)
attach(resList)

####################################################################################################
####################################################################################################
lbrks <- c(.0001,.0005,.001,.005,.01,.05,.1,.5,1)
cols <- punq$lab
names(cols) <- cols
cols <- c("NT"='red', "SWCT"='orange', "VR"='dark green', "RCT-gs"='dodger blue', "RCT"='dodger blue', "RCT-rp"='purple', "RCT-gs-rp" = 'purple')
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

## histogram of risk
jpeg(file.path(figdir, paste0('irsk.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(irsk[lab=='NT'], aes(x=inf)) + geom_histogram() + xlab('cumulative risk of infection') 
graphics.off()

## density lines risk spent
tmp <- irsk[type=='marg' & !lab%in%c('NT','VR')]
jpeg(file.path(figdir, paste0('irsk spentEV.jpeg')), w = wid, h = heig, units = 'in', res = res)
adj <- 2
p <- ggplot() + 
    geom_line(data=tmp, aes(x=spent_EV, col=lab, linetype = gs), stat='density', adjust=adj) +
              scale_color_manual(values=cols)  + xlab('per capita infection risk spent')  + ggtitle('average risk spent per individual') #+ xlim(.005, .3) + ylim(0,32)
print(p) ## print(p+scale_x_log10(breaks=lbrks))
graphics.off()

## infection risk
jpeg(file.path(figdir, paste0('irsk inf bars.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(irsk[lab=='NT'], aes(x=ordShow, y=inf, col=factor(cluster))) +  ylab('cumulative infection risk \n(over year after trial start)') + xlab('individual') +
    geom_bar(stat='identity', width=1) + ylim(0,1) + theme(legend.key.size = unit(.1, "cm")) + ggtitle('infection risk without vaccination') +
        theme(legend.position="top")
graphics.off()

## Conditional on arms & order randomization
tmp <- irsk[(arm==armShown & type=='cond' &  grepl('RCT',lab)) | (type=='condvd' & lab=='SWCT' & exmpl==T)]
####################################################################################################
## spent
ylim <- tmp[,range(spent,spent_EV)]
jpeg(file.path(figdir, paste0('irsk spent bars.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(tmp, aes(x=ordShowArm, y=spent, fill=armShown)) + ggtitle('risk spent, conditional on arm without EV') +
    geom_bar(stat='identity', width=1) + facet_wrap(~lab, ncol=2) +  ylab('risk') + xlab('individual') + ylim(ylim[1],ylim[2])
graphics.off()

jpeg(file.path(figdir, paste0('irsk spentEV bars.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(tmp, aes(x=ordShowArm, y=spent_EV, fill=armShown)) + ggtitle('risk spent, conditional on arm with EV') +
    geom_bar(stat='identity', width=1) + facet_wrap(~lab, ncol=2) +  ylab('risk')+ xlab('individual') + ylim(ylim[1],ylim[2])
graphics.off()

## avert
ylim <- tmp[,range(avert,avert_EV)]
jpeg(file.path(figdir, paste0('irsk avert bars.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(tmp, aes(x=ordShowArm, y=avert, fill=armShown)) + ggtitle('risk averted, conditional on arm without EV') +
    geom_bar(stat='identity', width=1) + facet_wrap(~lab, ncol=2) +  ylab('risk')+ xlab('individual') + ylim(ylim[1],ylim[2])
graphics.off()

jpeg(file.path(figdir, paste0('irsk avertEV bars.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(tmp, aes(x=ordShowArm, y=avert_EV, fill=armShown)) + ggtitle('risk averted, conditional on arm with EV') +
    geom_bar(stat='identity', width=1) + facet_wrap(~lab, ncol=2) +  ylab('risk')+ xlab('individual') + ylim(ylim[1],ylim[2])
graphics.off()
 
## both
ylim <- tmp[,range(-avert_EV,spent_EV, -avert, spent)]
jpeg(file.path(figdir, paste0('irsk spent & avert bars.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(tmp) + ggtitle('risk spent, conditional on arm without EV') +
    geom_hline(yintercept=.05) +
        geom_bar(aes(x=ordShowArm, y=spent, fill=armShown), stat='identity', width=1) +
            geom_bar(aes(x=ordShowArm, y=-avert, fill=armShown), stat='identity', width=1, alpha = .8) +
                facet_wrap(~lab, ncol=2) +  ylab('risk') + xlab('individual') + ylim(ylim[1],ylim[2])
graphics.off()

jpeg(file.path(figdir, paste0('irsk spent & avert EV bars.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(tmp) + ggtitle('risk spent, conditional on arm without EV') +
    geom_hline(yintercept=.05) +
        geom_bar(aes(x=ordShowArm, y=spent_EV, fill=armShown), stat='identity', width=1) +
            geom_bar(aes(x=ordShowArm, y=-avert_EV, fill=armShown), stat='identity', width=1, alpha = .8) +
                facet_wrap(~lab, ncol=2) +  ylab('risk') + xlab('individual') + ylim(ylim[1],ylim[2])
graphics.off()

####################################################################################################

## marginal on arms & order randomization
cols <- c('w/o EV' = 'dodger blue', 'w/ EV' = 'purple')
tmp <- irsk[type=='marg' & !lab %in% c('NT','VR')]
####################################################################################################
## spent
ylim <- tmp[,range(spent,spent_EV)]
jpeg(file.path(figdir, paste0('irsk spent bars MARG.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(tmp) + ggtitle('risk spent, marginal on randomization') +
        geom_bar(aes(x=ordShow, y=spent, fill=ev), data.table(tmp, ev = 'w/o EV'), stat='identity', width=1) +
        geom_bar(aes(x=ordShow, y=spent_EV, fill=ev), data.table(tmp[lab!='SWCT'], ev = 'w/ EV'), stat='identity', width=1) +
        scale_fill_manual(values=cols) +
        facet_wrap(~lab, ncol=2) +  ylab('risk') + xlab('individual') + ylim(ylim[1],ylim[2])
graphics.off()

## avert
ylim <- c(0, tmp[,max(avert,avert_EV)])
jpeg(file.path(figdir, paste0('irsk avert bars MARG.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(tmp) + ggtitle('risk averted, marginal on randomization') +
    geom_bar(aes(x=ordShow, y=avert_EV, fill=ev), data.table(tmp[lab!='SWCT'], ev = 'w/ EV'), stat='identity', width=1) +
    geom_bar(aes(x=ordShow, y=avert, fill=ev), data.table(tmp, ev = 'w/o EV'), stat='identity', width=1) +
        scale_fill_manual(values=cols) +
            facet_wrap(~lab, ncol=2) +
                ylab('risk')+ xlab('individual') + ylim(ylim[1],ylim[2])
graphics.off()

jpeg(file.path(figdir, paste0('haz traj.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(SpopH[simNum==1 & nbatch==3073], aes(Date, clusHaz*10^5, group = OcOrd, col=OcOrd)) + geom_line() + theme(legend.position="top") +
     theme(legend.key.size = unit(.1, "cm")) + ylab('daily infection hazard (per 100,000)') + ggtitle('mean cluster hazard trends')
graphics.off()

## **check why # of SWCT with posv is not the same as other designs** could be just running different sims?

jpeg(file.path(figdir, paste0('pow trhes.jpeg')), w = wid, h = heig, units = 'in', res = res)
ggplot(triSumm, aes(above_EV, power, col = lab)) + geom_point() + ylim(0,1) + xlab(paste0('expected # subjects spending >', threshold, ' infection risk'))
graphics.off()


## add error bars to inf spent & frac of information
## try fraction of information per person on y-axis
## show inf spent/averted on same plot

## equipoise perturbed plot
## histogram within strata groups
## look at risk by person/strata by treatment assignment for trials
## vaccinating a greater % of people increases risk averted/spent, but still has problem of withholding treatment from individuals
