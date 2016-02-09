if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('login1', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/work/EbolaVaccSim/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T)); gc()
require(optiRum); require(RColorBrewer); require(boot); require(data.table); require(ggplot2); require(grid); require(reshape2); require(parallel)
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('multiplot.R','extractFXN.R','ggplotTheme.R'), source)

thing <- 'Equip-indivL'
## Load VaccProp & hazT
load(file=paste0('BigResults/',thing,'/hazT',7,'.Rdata'))
load('data/vaccProp1.Rdata')
vaccProp <- vaccProp1
vaccProp[, simNum:=1:length(vaccEff)]
## Make Figure Folder
figdir <- file.path('Figures', thing)
dir.create(figdir)

## resList <- extractSims(thing, verb=0, maxbatches=NA, indivLev = T, mc.cores=48)
## resList <- procResList(resList)
## resList <- makeInfPow(resList)
load(file=file.path('BigResults',paste0(thing, '.Rdata')))
attach(resList)
names(resList)

ltys <- labsToShow <- finit[, levels(lab)]
names(ltys) <- ltys
ltys[c('RCT-gs-rp','SWCT')] <- 1; ltys[c('VR','NT')] <- 2; ltys[c('RCT','RCT-gs','RCT-rp')] <- 3 
class(ltys) <- 'numeric'
cols <- c("NT"='red', "RCT"='blue', "SWCT"='orange', "VR"='dark green', "RCT-gs"='magenta', "RCT-rp"='purple', "RCT-gs-rp" = 'dodger blue')

####################################################################################################
## set up hazard trajectories
rect <- data.table(xmin=as.Date(finit[,unique(trialStartDate)]), ymin=-Inf, ymax=Inf)
rect[,xmax :=xmin+24*7]
eb <- theme(
    ## axis.line=element_blank(),axis.text.x=element_blank(),
    #axis.text.y=element_blank(),
    axis.ticks=element_blank(),

    axis.title.x=element_blank(),
    axis.title.y=element_blank(),legend.position="none",
    panel.background=element_blank(),
    panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),plot.background=element_blank())
## HazT Plot
ip0 <- ggplot(hazT, aes(x=Date, y=clusHaz, col=as.factor(cluster))) + geom_line() + eb + scale_x_date(limits=as.Date(c('2014-08-01','2015-10-01'))) + labs(ylab='relative hazard')
ip <- ip0 + eb + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey20", border=NA, alpha=0.2, inherit.aes = FALSE)

####################################################################################################
## power vs infections averted

infAvertPow <- finit[grepl('Final',cat) & trial!='NT', list(.N, pow = mean(vaccGood), infAvert = mean(infAvert), 
                         infAvertProp = mean(infAvertProp)), ## infAvertableProp = mean(infAvertableProp, na.rm=T)), 
                     list(propInTrial, trialStartDate, cat, lab, avHaz)] 
infAvertPow[lab=='VR', pow:=0]
infAvertLab <- "infections averted relative to not performing any trial"

ftmp <- finit[grepl('Final',cat) & avHaz=='' & trial!='NT']
ftmp$catn <- factor(ftmp$cat)
ftmp[, catn:=factor(catn, labels = c('no rollout','rollout'))]
ftmp[,.N, list(propInTrial, trialStartDate, cat, lab, avHaz)]
tmp <- ftmp[lab %in% c('VR','NT')]
tmp[,cat:='allFinal_noEV']
tmp[,catn:='no rollout']
ftmp <- rbind(ftmp,tmp)

## density of infections averted by endtime, trial start date & % in trial
nmtmp <- 'infAvert dens .1 2014-10 Pos.pdf'
pdf(file.path(figdir, nmtmp), w=10, h = 8)
p <- ggplot(ftmp[posv==T & trialStartDate=='2014-10-01' & propInTrial==.1], aes(infAvert, colour = lab, linetype = lab)) +
    geom_density() + facet_grid(~catn) + 
        scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) + xlab(infAvertLab)
print(p)
graphics.off()

## density of infections averted by endtime, trial start date & % in trial
for(ii in 1:2) {
    if(ii==1) {
        ft <- ftmp[posv==T]
        nmtmp <- 'infAvert dens Pos.pdf'
    }else{
        ft <- ftmp
        nmtmp <- 'infAvert dens.pdf'
    }        
    pdf(file.path(figdir, nmtmp), w=10, h = 8)
    for(jj in 1:4) {
        ts <- ft[,unique(trialStartDate)][jj]
        p <- ggplot(ft[trialStartDate==ts], aes(infAvert, colour = lab, linetype = lab)) +
            geom_density() + labs(title=paste0('trial starts ', ts)) +
                ## facet_wrap(~propInTrial, ncol=1) +
                facet_grid(catn~propInTrial) +                 
                    scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) +
                                        #                scale_x_continuous(limits=c(-100,600)) + 
                        xlab(infAvertLab)
        ip1 <- ip + geom_rect(data=rect[jj], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey20", border=NA, alpha=0.2, inherit.aes = FALSE)
        multiplot(ip1, p, layout = matrix(c(1, rep(2,3)), ncol=1))
    }
    graphics.off()
}

ymaxHaz <- hazT[,max(clusHaz)]
rect[,ymax:=ymaxHaz*c(1,.9,.8,.7)]
ip1 <- ip0 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey20", border=NA, alpha=0.2, inherit.aes = FALSE)

## ethicline contours
tmp <- infAvertPow[cat=='allFinalEV' & lab!='NT' & avHaz=='']
ias <- pretty(tmp$infAvert_EV, n = 50)
pows <- seq(0,1,l=50)
contt <- data.table(expand.grid(ia=ias, pow=pows))
for(ii in 1:3) {
    infPpow <- c(100, 200, 300)[ii]
    contt[, merit:= ia/infPpow + pow] ## every infPpow infections is worth 10% power
    contt[pow<.3, merit:= ia/infPpow+.3] ## don't penalize less power under 30% power since it's so insignificant anyways
    ## Basic plot
    jpeg(file.path(figdir, paste0('infAvert pow', infPpow/10, '.jpeg')), w = 10, h = 8, units = 'in', res = 200)
    p <- ggplot() +
        ## fill
        geom_tile(data = contt, aes(x=ia, y=pow, z=merit, fill=merit)) + scale_fill_gradient(low = "beige", high = "brown") +
            stat_contour(data = contt, aes(x=ia, y=pow, z=merit), col = gray(.9, alpha = .9), bins = 7) +
                ## points
                geom_point(data = tmp, aes(infAvert_EV, pow, shape = lab, colour = lab, linetype = lab, size = 1.5)) + 
                    facet_grid(propInTrial~trialStartDate) + 
                        scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) + 
                            xlab(infAvertLab) + theme(legend.position='top', legend.box='horizontal') +
                                ## isoclines
                                geom_segment(data = tmp[lab=='VR'],
                                             aes(x = infAvert, y = .3, xend = 0, yend = .3+1/infPpow*infAvert)) +
                                                 geom_segment(data = tmp[lab=='VR'],
                                                              aes(x = infAvert, y = .3, xend = infAvert, yend = 0)) +  
                                                                  coord_cartesian(ylim=c(0,1)) + guides(size=F) +
                                                                      labs(title=paste(infPpow/10, 'infections := 10% power \n<30% power := negligible'))
    print(multiplot(ip1, p, layout = matrix(c(1, rep(2,3)), ncol=1)))
    graphics.off()
}

tmp <- infAvertPow[cat=='allFinalEV' & lab!='NT' & avHaz=='']
## by Proportion of infections averted
pdf(file.path(figdir, 'infAvertProp pow.pdf'), w = 10, h = 8)
    p <- ggplot(tmp, aes(infAvertProp, pow, shape = lab, colour = lab, linetype = lab)) +
        geom_point() + facet_grid(propInTrial~trialStartDate) + 
            scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) +
                        xlab(paste('proportion of ', infAvertLab)) +
                            geom_vline(aes(xintercept=infAvertProp), data = tmp[lab=='VR'], col = 'dark green')
multiplot(ip1, p, layout = matrix(c(1, rep(2,3)), ncol=1))
dev.off()

####################################################################################################
## individual level stuff

####################################################################################################
## fill in all hazard levels just so there aren't any missing ones, which can mess with things later
infAvertLab <- "infections averted relative to not performing any trial"

pdf(file.path(figdir, 'risk-stratified EV all trials & dates.pdf'), w = 10, h = 8)
tmp <- infAvertPow[avHaz=='' & trialStartDate=='2014-10-01'& !ihaz0Cat %in% levels(ihaz0Cat)[c(1:6,12:13)]]
p <- ggplot(tmp) + geom_point(aes(infSpent_EV, powfrac_EV, shape = lab, colour=lab, linetype=lab)) + facet_grid(ihaz0Cat~propInTrial) + 
    xlab('infections spent (EV)') + ylab('power fraction within strata')
print(p)
graphics.off()

####################################################################################################
## Compare RCT-gs-rp vs SWCT 
pdf(file.path(figdir, 'PoP vs infSpent RCT-gs-rp vs SWCT.pdf'), w = 8, h = 5)
tmp <- infAvertPow[avHaz=='' & !ihaz0Cat %in% levels(ihaz0Cat)[c(1:4,12:13)] & propInTrial==.05 & lab %in% c('SWCT','RCT-gs-rp')]
p <- ggplot(tmp) + geom_point(aes(infSpent_EV, powfrac_EV, colour=ihaz0Cat, shape=lab, size = 2)) + facet_grid(lab~trialStartDate) + 
    xlab('infections spent (EV)') + ylab('power fraction within strata') + scale_colour_brewer(palette='Reds') 
print(p)
graphics.off()

## no end vaccination
pdf(file.path(figdir, 'PoP_noEV vs infSpent_noEV RCT-gs-rp vs SWCT.pdf'), w = 8, h = 5)
p <- ggplot(tmp) + geom_point(aes(infSpent_noEV, powfrac_noEV, colour=ihaz0Cat, shape=lab, size = 2)) + facet_wrap(~trialStartDate, nc=4) + 
    xlab('infections spent (EV)') + ylab('power fraction within strata') + scale_colour_brewer(palette='Reds') 
print(p)
graphics.off()
####################################################################################################

infAvertPow[,hazLab:=as.numeric(ihaz0Cat)]
tmp0 <- infAvertPow[avHaz=='' & propInTrial==.05 & lab %in%  c('SWCT','RCT-gs-rp')]
tmp <- rbind(tmp0[, list(trialStartDate, hazLab, lab, avHaz, type = 'avert', val = infAvertPC_EV)],
             tmp0[, list(trialStartDate, hazLab, lab, avHaz, type = 'spent', val = infSpentPC_EV)])
pdf(file.path(figdir, 'per capita risk averted by strata.pdf'), w = 8, h = 5)
p <- ggplot(tmp) + geom_bar(aes(hazLab, val, fill = type), stat='identity', position='dodge') + facet_wrap(lab~trialStartDate, nc=4)  +
    xlab('per capita risk') + xlab('hazard level') + scale_color_manual(values=c('red','blue'))
print(p)
graphics.off()
## ** CHANGED extractFXn to haz0cat, should change back to ihaz0cat ****


jpeg(file.path(figdir, paste0('infSpentPC_EV versus PoP.jpeg')), w = 10, h = 8, units = 'in', res = 200)
tmp <- infAvertPow[avHaz=='' & propInTrial==.1 & !ihaz0Cat %in% levels(ihaz0Cat)[c(1:4,13)] & lab!='VR']
p <- ggplot() +
    geom_point(data = tmp, aes(infSpentPC_EV, powfrac_EV, shape = lab, colour = ihaz0Cat, size = 1.5)) + 
        facet_grid(~trialStartDate) + scale_colour_manual(values=rev(brewer.pal(n=length(unique(tmp$ihaz0Cat)), "Spectral"))) +
                xlab('per capita risk spent') + ylab('fraction of information from hazard class') 
print(p)
graphics.off()

## proportion of avertable infections averted
jpeg(file.path(figdir, paste0('proportion avertable averted versus PoP.jpeg')), w = 10, h = 8, units = 'in', res = 200)
tmp <- infAvertPow[avHaz=='' & propInTrial==.1 & !ihaz0Cat %in% levels(ihaz0Cat)[c(1:4,13)] & lab!='VR']
tmp[,pAvert:=infAvert_EV/(infAvert_EV+infSpent_EV)]
p <- ggplot() +
    geom_point(data = tmp, aes(pAvert, powfrac_EV, shape = lab, colour = ihaz0Cat, size = 1.5)) + 
        facet_grid(~trialStartDate) + scale_colour_manual(values=rev(brewer.pal(n=length(unique(tmp$ihaz0Cat)), "Spectral"))) +
                xlab('proportion avertable infections averted') + ylab('fraction of information from hazard class')  + xlim(0,1)
print(p)
graphics.off()

jpeg(file.path(figdir, paste0('infSpent_EV pow.jpeg')), w = 10, h = 8, units = 'in', res = 200)
tmp <- infAvertPow[avHaz=='' & propInTrial==.1 & !ihaz0Cat %in% levels(ihaz0Cat)[c(1:2,12:13)]]
p <- ggplot() +
    geom_point(data = tmp, aes(pow_per_infSpent_EV, power, shape = lab, colour = lab, linetype = lab, size = 1.5)) + 
        facet_grid(ihaz0Cat~trialStartDate) + 
            scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) + xlim(-.1,.4) + 
                xlab('power / infection not averted') + theme(legend.position='top', legend.box='horizontal') ## +
print(p)
graphics.off()


tmp <- infAvertPow[avHaz=='' & propInTrial==.1 & !ihaz0Cat %in% levels(ihaz0Cat)[c(1:4,13)]]
## Basic plot
jpeg(file.path(figdir, paste0('infAvertPC_EV pow.jpeg')), w = 10, h = 8, units = 'in', res = 200)
p <- ggplot() +
    ## points
    geom_point(data = tmp, aes(infAvertPC_EV, powfrac_EV, shape = lab, colour = lab, linetype = lab, size = 1.5)) + 
        facet_grid(ihaz0Cat~trialStartDate) + 
            scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) + xlim(-.05,.15) + 
                xlab('per capita risk averted') + ylab('fraction of information from hazard class') + theme(legend.position='top', legend.box='horizontal') ## +
print(p)

graphics.off()

## per capita incidence averted & spent always add up to same value for given strata across designs.
rowSums(tmp[trialStartDate=='2014-10-01' & ihaz0Cat=='(0.01,0.1]', list(infAvertPC_noEV, infSpentPC_noEV)])

## add error bars to inf spent & frac of information
## try fraction of information per person on y-axis
## show inf spent/averted on same plot
