if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T)); gc()
require(RColorBrewer); require(boot); require(data.table); require(ggplot2); require(grid); require(reshape2); require(parallel)
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('multiplot.R','extractFXN.R','ggplotTheme.R'), source)

####################################################################################################
## extract factuals
thing <- 'Equip-indivL'
out <- extractSims(thing, verb=0, maxbatches=NA, indivLev = T, mc.cores=48)
load(file=file.path('BigResults',paste0(thing, '.Rdata')))
finTrials[order(gs), list(tcalMean = mean(tcal), power = mean(vaccGood), length(tcal)),
          list(trial, gs, ord, delayUnit, propInTrial, trialStartDate, avHaz)]


## finInfo has 4 categories for each simulation in finTrials (all=all cases by trial end date,
## analyzed = all analyzed cases by trial end date, allFinalEV/no_EV = all cases by stop date w/ or
## w/o end vaccine rollout)
load(file=paste0('BigResults/',thing,'/hazT',7,'.Rdata'))

figdir <- file.path('Figures', thing)
dir.create(figdir)

####################################################################################################
## extract counterfactuals
load('data/vaccProp1.Rdata')
vaccProp <- vaccProp1
vaccProp[, simNum:=1:length(vaccEff)]


setkey(finInfo, gs, ord, trial, propInTrial, simNum, trialStartDate, cat)
setkey(finTrials, gs, ord, trial, propInTrial, simNum, trialStartDate)

## Merge them so we have simulation level-info (speed, result, etc) in each finInfo category
finit <- finInfo[finTrials[, list(gs, ord, trial, propInTrial, trialStartDate, simNum, tcal, vaccCases, contCases, vaccGood, vaccBad, vaccEff, pSAE, PHU)]]
finit[, lab:=factor(paste0(trial, c('','gs-')[as.numeric(gs==T)+1],'-', ord))] ## make useful labels

## Merge finit with fincfs so we can compare counterfactuals and factuals
fall <- rbind(fincfs[, list(lab = cf, cat = 'allFinalEV', propInTrial, caseTot, vaccEff, simNum, trialStartDate, nsae, pSAE)],
              fincfs[, list(lab = cf, cat = 'allFinal_noEV', propInTrial, caseTot, vaccEff, simNum, trialStartDate, nsae, pSAE)],
              finit[, list(lab, cat, propInTrial, caseTot, vaccEff, simNum, trialStartDate, nsae, pSAE, vaccGood)], fill=T) ##cat %in% c('all','allFinalEV','allFinal_noEV')
fall[, posv:= vaccEff>0] ## positive vaccine efficacy simulations

fall[, list(nsim = length(caseTot), meansae = mean(nsae)), list(lab, cat, propInTrial, trialStartDate)]##[,range(V1)]

fall[, unique(lab)]
fall[grepl('Final',cat), caseTotNT := caseTot[lab=='NTpop'], list(simNum, propInTrial, trialStartDate)]
fall[, infAvert := caseTotNT - caseTot]
fall[, infAvertprop := infAvert/caseTotNT]
fall[, infAvertSAE := caseTotNT - (caseTot+nsae)]
fall[, infAvertSAEprop := infAvertSAE/caseTotNT]

labsToShow <- c('RCTgs--TU','SWCT-none', 'VRpop','RCT-none')
ltys <- c('VRpop'=2, 'RCT-none' = 3, 'RCTgs--TU'=1, 'SWCT-none'=1, 'NTpop'=2)
cols <- c('VRpop'='dark green', 'RCT-none' = 'purple', 'RCTgs--TU'="#333BFF", 'SWCT-none'='orange', 'NTpop' = 'red')

ftmp <- fall[lab %in% labsToShow & grepl('Final',cat)]
ftmp$catn <- factor(ftmp$cat)
ftmp[, catn:=factor(catn, labels = c('1yr w/o rollout','1yr w/ rollout'))]

####################################################################################################
## set up hazard trajectories
rect <- data.table(xmin=as.Date(infAvertPow[,unique(trialStartDate)]), ymin=-Inf, ymax=Inf)
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
ip0 <- ggplot(hazT, aes(x=Date, y=clusHaz, col=as.factor(cluster))) + geom_line() + eb + scale_x_date(limits=as.Date(c('2014-08-01','2015-10-01'))) + labs(ylab='relative hazard')
####################################################################################################

####################################################################################################
## power vs infections averted
ip <- ip0 + eb + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey20", border=NA, alpha=0.2, inherit.aes = FALSE)

infAvertPow <- fall[!lab %in% c('NTpop'), list(pow = mean(vaccGood), infAvert = mean(infAvert), infAvertSAE = mean(infAvertSAE),
                                                  infAvertprop = mean(infAvertprop), infAvertSAEprop = mean(infAvertSAEprop)), list(propInTrial, trialStartDate, cat, lab)]
infAvertPow[lab=='VRpop', pow:=0]
infAvertPow[,infAvertpropVR := infAvertprop/infAvertprop[lab=='VRpop'], list(propInTrial,trialStartDate,cat)]

infAvertLab <- "infections averted relative to not performing any trial"

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
                                ##                        geom_vline(aes(xintercept=infAvert), data = infAvertPow[lab=='VRpop'], col = 'dark green') +
                                geom_segment(data = infAvertPow[lab=='VRpop'&cat=='allFinalEV'],
                                             aes(x = infAvert, y = .3, xend = 0, yend = .3+1/infPpow*infAvert)) +
                                                 geom_segment(data = infAvertPow[lab=='VRpop'&cat=='allFinalEV'],
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
                            geom_vline(aes(xintercept=infAvertprop), data = infAvertPow[lab=='VRpop'], col = 'dark green')
multiplot(ip1, p, layout = matrix(c(1, rep(2,3)), ncol=1))
dev.off()

## by proportion of avertable infections averted
pdf(file.path(figdir, 'infAvertpropVR pow.pdf'), w = 10, h = 8)
p <- ggplot(infAvertPow[cat=='allFinalEV' & lab %in% labsToShow], aes(infAvertpropVR, pow, colour = lab, linetype = lab)) +
    geom_point() + facet_grid(propInTrial~trialStartDate) + 
        scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) +
            xlab(paste('proportion of avertable (from VR) deaths averted')) +
                geom_vline(aes(xintercept=infAvertpropVR), data = infAvertPow[lab=='VRpop'], col = 'dark green')
multiplot(ip1, p, layout = matrix(c(1, rep(2,3)), ncol=1))
dev.off()



####################################################################################################
## Figures
####################################################################################################
labsToShow <- c('RCTgs--TU','SWCT-none', 'VRpop','RCT-none', 'NTpop')
ftmp <- fall[lab %in% labsToShow & grepl('Final',cat)]
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
