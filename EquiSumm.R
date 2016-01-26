if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T)); gc()
require(RColorBrewer); require(boot); require(data.table); require(ggplot2); require(grid); require(reshape2); require(parallel)
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
## load(file=file.path('BigResults',paste0(thing, '.Rdata')))
## attach(resList)

## parms[,c('nbsize','weeklyDecay','cvWeeklyDecay','cvClus','cvClusTime'):=NULL]
## finInfo[,1:9,with=F]
## finTrials[,1:9,with=F]

## ## merge all trial-level info
## finit <- merge(parms, merge(finInfo, finTrials, all.x=T), all.y=T, by = 'nbatch')

## finit[order(gs), list(tcalMean = mean(tcal), power = mean(vaccGood), length(tcal)),
##           list(trial, gs, ord, delayUnit, propInTrial, trialStartDate, avHaz, cat)]
## cols <- c('trial', 'gs', 'ord', 'delayUnit', 'propInTrial', 'trialStartDate', 'avHaz', 'cat')
## finit[, (cols):=lapply(.SD, as.factor), .SDcols=cols]
## ## make sure all simulations completed (2040 of each)
## nsms <- finit[, .N, list(trial, gs, ord, delayUnit, propInTrial, trialStartDate, avHaz, cat)]
## cols <- c('trial', 'gs', 'ord', 'delayUnit', 'propInTrial', 'trialStartDate', 'avHaz', 'cat')
## nsms[, (cols):=lapply(.SD, as.factor), .SDcols=cols]
## summary(nsms[N<2040]) ## SWCT's are missing results, working on this

## finit[, posv:= vaccEff>0] ## positive vaccine efficacy simulations
## names(finit)
## finit[,list(trial, sim, simNum, cat, caseTot)]

## ## Figure out how to match trials
## finit[,.N,list(nbatch, sim, simNum, cat)] ## full unique
## finit[trial=='NT',.N,list(propInTrial, trialStartDate, avHaz, sim, simNum)] ## unique CF
## finit[trial=='VR',.N,list(propInTrial, trialStartDate, avHaz, sim, simNum)] ## unique CF
## finit[cat=='allFinalEV',.N,list(propInTrial, trialStartDate, avHaz, sim, simNum)] ## matched scenarios (i.e. 7 factuals/2cfs)

## finit[, c('caseTotNT','caseTotVR'):=as.numeric(NA)]
## ## EV
## finit[cat=='allFinalEV', caseTotNT:=caseTot[trial=='NT'], list(propInTrial, trialStartDate, avHaz, sim, simNum)]
## finit[cat=='allFinalEV', caseTotVR:=caseTot[trial=='VR'], list(propInTrial, trialStartDate, avHaz, sim, simNum)]
## ## noEV
## finit[(trial %in% c('RCT','SWCT') & cat=='allFinal_noEV')|(trial=='NT' & cat=='allFinalEV'), caseTotNT:=caseTot[trial=='NT'],
##       list(propInTrial, trialStartDate, avHaz, sim, simNum)] #     
## finit[(trial %in% c('RCT','SWCT') & cat=='allFinal_noEV')|(trial=='VR' & cat=='allFinalEV'), caseTotVR:=caseTot[trial=='VR'],
##       list(propInTrial, trialStartDate, avHaz, sim, simNum)] #     
## ## inf Avert
## finit[, infAvert := caseTotNT - caseTot]
## finit[, list(propInTrial, trialStartDate, avHaz, sim, simNum, trial, gs, ord, cat, caseTot, caseTotNT, caseTotVR, infAvert)] #     
## finit[, infAvertProp := infAvert/caseTotNT]
## finit[, infAvertableProp := infAvert/(caseTotNT-caseTotVR)]
## finit[, lab:=trial]
## finit[trial=='RCT' & gs==T, lab:=paste0(lab,'-gs')]
## finit[trial=='RCT' & ord=='TU', lab:=paste0(lab,'-rp')]
## finit[,unique(lab)]
## finit[,lab:=as.factor(lab)]
## save(finit, file=file.path('BigResults',paste0(thing, '_finit.Rdata')))
load(file=file.path('BigResults',paste0(thing, '_finit.Rdata')))

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
                         infAvertProp = mean(infAvertProp)), #infAvertableProp = mean(infAvertableProp)), 
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
    p <- ggplot(ft[trialStartDate=='2014-10-01' & propInTrial==.1], aes(infAvert, colour = lab, linetype = lab)) +
        geom_density() + facet_grid(~catn) + 
            scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) +
                                        #                scale_x_continuous(limits=c(-100,600)) + 
                xlab(infAvertLab)
    ip1 <- ip + geom_rect(data=rect[jj], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey20", border=NA, alpha=0.2, inherit.aes = FALSE)
    multiplot(ip1, p, layout = matrix(c(1, rep(2,3)), ncol=1))
    graphics.off()
}

finit[avHaz=='' & cat=='allFinal_noEV' & trial=='RCT' & trialStartDate=='2014-10-01' & propInTrial==.1 & sim ==85 & simNum==2040, 
      list(nbatch, sim, simNum, trial, ord, gs, infAvert, vaccEff) ]

parmsMat[c(792,816,840,864)]

finit[nbatch==792 & sim==85 & simNum==2040 & trial=='RCT']

ymaxHaz <- hazT[,max(clusHaz)]
rect[,ymax:=ymaxHaz*c(1,.9,.8,.7)]
ip1 <- ip0 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey20", border=NA, alpha=0.2, inherit.aes = FALSE)

## ethicline contours
tmp <- infAvertPow[cat=='allFinalEV' & lab %in% labsToShow]
ias <- pretty(tmp$infAvert, n = 50)
pows <- seq(0,1,l=50)
contt <- data.table(expand.grid(ia=ias, pow=pows))
for(ii in 1:3) {
    infPpow <- c(100, 200, 300)[ii]
    contt[, merit:= ia/infPpow + pow] ## every infPpow infections is worth 10% power
    contt[pow<.3, merit:= ia/infPpow+.3] ## don't penalize less power under 30% power since it's so insignificant anyways
    v <- ggplot(contt, aes(ia,pow,z=merit)) +
        geom_tile(aes(fill=merit))+ stat_contour() + scale_fill_gradient(low = "brown", high = "white")

    ## Basic plot
    jpeg(file.path(figdir, paste0('infAvert pow', infPpow/10, '.jpeg')), w = 10, h = 8, units = 'in', res = 200)
    p <- ggplot() +
        geom_tile(data = contt, aes(x=ia, y=pow, z=merit, fill=merit)) + scale_fill_gradient(low = "brown", high = "pink") +
            stat_contour(data = contt, aes(x=ia, y=pow, z=merit), col = gray(.9, alpha = .9), bins = 7) +
                geom_point(data = infAvertPow[cat=='allFinalEV' & lab %in% labsToShow], aes(infAvert, pow, colour = lab, linetype = lab, size = 1.5)) + 
                    facet_grid(propInTrial~trialStartDate) + 
                        scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) + 
                            xlab(infAvertLab) + theme(legend.position='top', legend.box='horizontal') +
                                ##                        geom_vline(aes(xintercept=infAvert), data = infAvertPow[lab=='VR'], col = 'dark green') +
                                geom_segment(data = infAvertPow[lab=='VR'&cat=='allFinalEV'],
                                             aes(x = infAvert, y = .3, xend = 0, yend = .3+1/infPpow*infAvert)) +
                                                 geom_segment(data = infAvertPow[lab=='VR'&cat=='allFinalEV'],
                                                              aes(x = infAvert, y = .3, xend = infAvert, yend = 0)) +  
                                                                  coord_cartesian(ylim=c(0,1)) + guides(size=F) +
                                                                      labs(title=paste(infPpow/10, 'infections := 10% power \n<30% power := negligible'))
    print(multiplot(ip1, p, layout = matrix(c(1, rep(2,3)), ncol=1)))
    graphics.off()
}

## by proportion of infections averted
pdf(file.path(figdir, 'infAvertprop pow.pdf'), w = 10, h = 8)
    p <- ggplot(infAvertPow[cat=='allFinalEV' & lab %in% labsToShow], aes(infAvertprop, pow, colour = lab, linetype = lab)) +
        geom_point() + facet_grid(propInTrial~trialStartDate) + 
            scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) +
                        xlab(paste('proportion of ', infAvertLab)) +
                            geom_vline(aes(xintercept=infAvertprop), data = infAvertPow[lab=='VR'], col = 'dark green')
multiplot(ip1, p, layout = matrix(c(1, rep(2,3)), ncol=1))
dev.off()

## by proportion of avertable infections averted
pdf(file.path(figdir, 'infAvertpropVR pow.pdf'), w = 10, h = 8)
p <- ggplot(infAvertPow[cat=='allFinalEV' & lab %in% labsToShow], aes(infAvertpropVR, pow, colour = lab, linetype = lab)) +
    geom_point() + facet_grid(propInTrial~trialStartDate) + 
        scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) +
            xlab(paste('proportion of avertable (from VR) deaths averted')) +
                geom_vline(aes(xintercept=infAvertpropVR), data = infAvertPow[lab=='VR'], col = 'dark green')
multiplot(ip1, p, layout = matrix(c(1, rep(2,3)), ncol=1))
dev.off()



####################################################################################################
## Figures
####################################################################################################
labsToShow <- c('RCTgs--TU','SWCT-none', 'VR','RCT-none', 'NT')
ftmp <- finit[lab %in% labsToShow & grepl('Final',cat)]
ftmp$catn <- factor(ftmp$cat)
ftmp[, catn:=factor(catn, labels = c('1yr w/o rollout','1yr w/ rollout'))]

ty <- theme(axis.text.y=element_blank(),
            axis.title.y=element_blank())

## pit facets, tsd pages
pdf(file.path(figdir, 'case dens.pdf'))
for(jj in 1:4) {
    ts <- ft[,unique(trialStartDate)][jj]
    p <- ggplot(ftmp[trialStartDate==ts], aes(caseTot, colour = lab, linetype = lab)) +
        geom_density() + facet_wrap(~propInTrial, ncol=1, scales='free_y') + 
            scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) +
            xlab("total infections by 1 year post start date")
    ip2 <- ip + geom_rect(data=rect[jj], aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey20", border=NA, alpha=0.2, inherit.aes = FALSE)
    multiplot(ip2, p+ty, layout = matrix(c(1, rep(2,3)), ncol=1))
}
dev.off()

ymaxHaz <- hazT[,max(clusHaz)]
rect[,ymax:=ymaxHaz*c(1,.9,.8,.7)]
ip1 <- ip + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey20", border=NA, alpha=0.2, inherit.aes = FALSE)
## pit pages, tsd facets
pdf(file.path(figdir, 'caseXpitXtsd.pdf'))
for(jj in 1:4) {
    ts <- ft[,unique(propInTrial)][jj]
    p <- ggplot(ftmp[propInTrial==ts], aes(caseTot, colour = lab, linetype = lab)) +
        geom_density() + facet_wrap(~trialStartDate, ncol=1, scales='free_y') + labs(title=paste0('proportion cases in trial ', ts)) +
            scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) +
                xlab("total infections by 1 year post start date")
    multiplot(ip1, p+ty, layout = matrix(c(1, rep(2,3)), ncol=1))
}
dev.off()


## break down what's happening before trial & after trial (due to roll out)
## parse out false pos/negative from correct identifications
## is there a simple linear relationship between equipoise & power? (medium risk vs high risk group, do you lose more info than equipoise gained)
## post-vacc rollout period may modulate linearity of this relationship
