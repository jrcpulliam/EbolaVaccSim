if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
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
ias <- pretty(tmp$infAvert, n = 50)
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
                geom_point(data = tmp, aes(infAvert, pow, shape = lab, colour = lab, linetype = lab, size = 1.5)) + 
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

shc <- CJ.dt(unique(strat_haz0Cat[,list(nbatch, simNum)]), data.table(haz0Cat = strat_haz0Cat[,unique(haz0Cat)]))
setkey(shc, nbatch, simNum, haz0Cat)
setkey(strat_haz0Cat, nbatch, simNum, haz0Cat)
setkey(shc, nbatch, simNum, haz0Cat)
shc <- merge(shc, strat_haz0Cat, all=T) ##
shc[is.na(N_noEVinf), N_noEVinf:=0] ## replace NA's with 0s
shc[is.na(N_noEVsae), N_noEVsae:=0]
shc[!nbatch %in% parms[trial %in% c('NT','VR'), nbatch] & is.na(N_EVinf), N_EVinf:=0] ## replace NAs w 0s only if a factual (otherwise should be NA
shc[!nbatch %in% parms[trial %in% c('NT','VR'), nbatch] & is.na(N_EVsae), N_EVsae:=0]
setkey(shc, nbatch, simNum, haz0Cat)

shc <- merge(parms, shc, all.y=T, by = c('nbatch'))
shc <- merge(shc, finTrials[,list(nbatch, simNum, vaccGood, tcal)], by = c('nbatch','simNum'))
setkey(shc, nbatch, simNum, haz0Cat)


shc[,.N, list(propInTrial, ord, gs, trialStartDate, trial, simNum, nbatch, avHaz)] # number of hazard categories
shc[,.N, list(propInTrial, ord, gs, trialStartDate, trial, simNum, nbatch, avHaz, haz0Cat)] # down to unique

shc[, .N, list(propInTrial, trialStartDate, simNum, avHaz, haz0Cat)] ## 7 simulation types
unique(parms[propInTrial==.025 & trialStartDate=='2014-10-01' & avHaz==''][,2:12,with=F])

shc[, N_NT:= N_noEVinf[trial=='NT'], list(propInTrial, trialStartDate, simNum, avHaz, haz0Cat)]
shc[, N_VR:= N_noEVinf[trial=='VR'], list(propInTrial, trialStartDate, simNum, avHaz, haz0Cat)]
shc[, infAvert := N_NT - N_EVinf]
shc[, infAvert_noEV := N_NT - N_noEVinf]
shc[, infAvertProp := infAvert/N_NT]
shc[, infAvert_noEVProp := infAvert_noEV/N_NT]

shc[, lab:=trial]
shc[trial=='RCT' & gs==T, lab:=paste0(lab,'-gs')]
shc[trial=='RCT' & ord=='TU', lab:=paste0(lab,'-rp')]
shc[,lab:=as.factor(lab)]

infAvertPow <- shc[trial!='NT', list(.N, infAvert = mean(infAvert), infAvertProp = mean(infAvertProp)
                   , infAvert_noEV = mean(infAvert_noEV), infAvert_noEVProp = mean(infAvert_noEVProp)), 
                     list(propInTrial, trialStartDate, lab, avHaz, haz0Cat)]

infAvertPow <- merge(infAvertPow
    , finit[cat=='allFinalEV' & trial!='NT', list(.N, power=mean(vaccGood)), list(propInTrial, trialStartDate, lab, avHaz)]
    , by = c('propInTrial', 'trialStartDate', 'lab', 'avHaz'))
infAvertPow[lab=='VR', power:=0]

infAvertLab <- "infections averted relative to not performing any trial"

pdf(file.path(figdir, 'riskstrat.pdf'), w = 10, h = 8)
p <- ggplot(shc[avHaz=='' & trialStartDate=='2014-10-01']) + 
geom_point(aes(N_EVinf

graphics.off()



