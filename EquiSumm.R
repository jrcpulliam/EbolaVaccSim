if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T))
require(RColorBrewer); require(boot); require(data.table); require(vioplot); require(ggplot2); require(grid); require(reshape2)
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('simFuns.R','AnalysisFuns.R','CoxFxns.R','EndTrialFuns.R', 'extractFXN.R'), source)
source('ggplotTheme.R')
source('data/plotIncCountry.R')
source('multiplot.R')

####################################################################################################
## extract factuals
thing <- 'Equip-ByTrialDate' ## 'Equip-rand3pit'
## out <- extractSims(thing, verb=0)
load(file=file.path('BigResults',paste0(thing, '.Rdata')))
finTrials[order(gs), list(tcalMean = mean(tcal), power = mean(vaccGood), length(tcal)),
          list(trial, gs, ord, delayUnit, propInTrial, trialStartDate)]
finTrials[tcal<168 & vaccEff>0, list(nsim=length(tcal), tcalMean = mean(tcal), power = mean(vaccGood)),
          list(trial, gs, ord, delayUnit)]
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

thingCFs <- paste0(thing, 'CFs') ##'Equip-cfs-3pit'
## fincfs <- extractCFs(thingCFs, verb=0)
load(file=file.path('BigResults',paste0(thingCFs, '.Rdata')))
fincfs <- merge(fincfs, vaccProp, by = 'simNum', all.y=F) ## copy vaccProp in there (should do this in analysisFuns.R later
class(fincfs$cc) <- 'numeric'
setnames(fincfs, 'saeTot','nsae')

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
ltys <- c('VRpop'=2, 'RCT-none' = 3, 'RCTgs--TU'=1, 'SWCT-none'=1)
cols <- c('VRpop'='dark green', 'RCT-none' = 'purple', 'RCTgs--TU'="#333BFF", 'SWCT-none'='orange')

ftmp <- fall[lab %in% labsToShow & grepl('Final',cat)]
ftmp$catn <- factor(ftmp$cat)
ftmp[, catn:=factor(catn, labels = c('1yr w/o rollout','1yr w/ rollout'))]

####################################################################################################
## power vs infections averted
ip <- ggplot(hazT, aes(x=Date, y=clusHaz, col=as.factor(cluster))) + geom_line() + eb + scale_x_date() +
    geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey20", border=NA, alpha=0.2, inherit.aes = FALSE)

infAvertPow <- fall[!lab %in% c('NTpop'), list(pow = mean(vaccGood), infAvert = mean(infAvert), infAvertSAE = mean(infAvertSAE),
                                                  infAvertprop = mean(infAvertprop), infAvertSAEprop = mean(infAvertSAEprop)), list(propInTrial, trialStartDate, cat, lab)]
infAvertPow[lab=='VRpop', pow:=0]
infAvertPow[,infAvertpropVR := infAvertprop/infAvertprop[lab=='VRpop'], list(propInTrial,trialStartDate,cat)]

infAvertLab <- "infections averted relative to not performing any trial"

####################################################################################################
## set up hazard trajectories
ymax <- hazT[,max(clusHaz)]
rect <- data.table(xmin=as.Date(infAvertPow[,unique(trialStartDate)]), ymin=-Inf, ymax=Inf)
rect[,xmax :=xmin+24*7]
eb <- theme(
    ## axis.line=element_blank(),axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),legend.position="none",
    panel.background=element_blank(),
    panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),plot.background=element_blank())
ip <- ggplot(hazT, aes(x=Date, y=clusHaz, col=as.factor(cluster))) + geom_line() + eb + scale_x_date() ##+
####################################################################################################

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

rect[,ymax:=ymax*c(1,.9,.8,.7)]
ip <- ip + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey20", border=NA, alpha=0.2, inherit.aes = FALSE)

## Basic plot
pdf(file.path(figdir, 'infAvert pow.pdf'), w = 10, h = 8)
p <- ggplot(infAvertPow[cat=='allFinalEV' & lab %in% labsToShow], aes(infAvert, pow, colour = lab, linetype = lab)) +
    geom_point() + facet_grid(propInTrial~trialStartDate) + 
        scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) +
            xlab(infAvertLab) + theme(legend.position='top')
                geom_vline(aes(xintercept=infAvert), data = infAvertPow[lab=='VRpop'], col = 'dark green')
multiplot(ip, p, layout = matrix(c(1, rep(2,3)), ncol=1))
dev.off()

## by proportion of infections averted
pdf(file.path(figdir, 'infAvertprop pow.pdf'), w = 10, h = 8)
    p <- ggplot(infAvertPow[cat=='allFinalEV' & lab %in% labsToShow], aes(infAvertprop, pow, colour = lab, linetype = lab)) +
        geom_point() + facet_grid(propInTrial~trialStartDate) + 
            scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) +
                        xlab(paste('proportion of ', infAvertLab)) +
                            geom_vline(aes(xintercept=infAvertprop), data = infAvertPow[lab=='VRpop'], col = 'dark green')
multiplot(ip, p, layout = matrix(c(1, rep(2,3)), ncol=1))
dev.off()

## by proportion of avertable infections averted
pdf(file.path(figdir, 'infAvertpropVR pow.pdf'), w = 10, h = 8)
p <- ggplot(infAvertPow[cat=='allFinalEV' & lab %in% labsToShow], aes(infAvertpropVR, pow, colour = lab, linetype = lab)) +
    geom_point() + facet_grid(propInTrial~trialStartDate) + 
        scale_color_manual(values=cols) + scale_linetype_manual(values = ltys) +
            xlab(paste('proportion of avertable (from VR) deaths averted'))
multiplot(ip, p, layout = matrix(c(1, rep(2,3)), ncol=1))
dev.off()
              
####################################################################################################
## Figures
####################################################################################################
fall[,unique(lab)]


## pit facets, tsd pages
pdf(file.path(figdir, 'case dens.pdf'))
for(ts in ftmp[,unique(trialStartDate)]) {
    p <- ggplot(ftmp[trialStartDate==ts], aes(caseTot, colour = lab, linetype = lab)) +
        geom_density() + facet_wrap(~propInTrial, ncol=1) + labs(title=paste0('trial starts ', ts)) +
            scale_color_manual(values=groupcols) + scale_linetype_manual(values = ltys) +
                xlab("total cases by 1 year post start date")
    print(p+thsb)
}
dev.off()

## pit pages, tsd facets
pdf(file.path(figdir, 'caseXpitXtsd.pdf'))
for(ts in ftmp[,unique(propInTrial)]) {
    p <- ggplot(ftmp[propInTrial==ts], aes(caseTot, colour = lab, linetype = lab)) +
        geom_density() + facet_wrap(~trialStartDate, ncol=1) + labs(title=paste0('proportion cases in trial ', ts)) +
            scale_color_manual(values=groupcols) + scale_linetype_manual(values = ltys) +
                xlab("total cases by 1 year post start date")
    print(p+thsb)
}
dev.off()



####################################################################################################
## Difference between matched counterfactuals and factuals

## Merge (big merger but worth it in speed later)
fin <- merge(finit, fincfs, by ='simNum', all.y=F, allow.cartesian = T)
class(fin$simNum) <- class(fin$cc)<- 'integer'
setkey(fin, gs, ord.x, simNum, cf, cat)

## Get difference between caseTot
ctot <- fin[cat %in% c('allFinalEV', 'allFinal_noEV'), 
            list(cc, vaccGood, vaccBad, tcal, caseTotvsCF = caseTot.x - caseTot.y, caseTotF = caseTot.x, caseTotCF = caseTot.y),             list(gs, ord.x, simNum, cf, cat)]
setkey(ctot, simNum, cc, cf)
satr <- CJ(simNum = 1:max(ctot$simNum), cc = 1:max(ctot$cc), cf = c('NTpop','VRpop'))
ctot <- ctot[satr, allow.cartesian=T]## Merge to make sure that empty columns exist
ctot

## SEGAULT HERE
## Get difference between caseTot
## [1] below is because rows are replicated by cf, except for cf & caseTot.y & SAE
ctot <- fin[cat %in% c('allFinalEV', 'allFinal_noEV'), 
            list(cc[1], vaccGood[1], vaccBad[1], tcal[1],
                 caseTotF = caseTot.x[1],
                 caseTotNT = caseTot.y[cf=='NTpop'],
                 caseTotVR = caseTot.y[cf=='VRpop']), 
            list(gs, ord.x, simNum, cat)]
ctot

## ctot[, list(ctF = mean(caseTotF),
##             ctCF = mean(caseTotCF),
##             diff = mean(caseTotvsCF)),
##      list(gs, ord.x, simNum, cf, cat)]

## cs <- ctot[order(cat), list(ctF = mean(caseTotF), ctCF = mean(caseTotCF), diff = mean(caseTotvsCF), 
##                             power = mean(vaccGood), tcal = mean(tcal)), 
##            list(gs, ord.x, cf, cat)]
## save(cs, file = file.path('Results','cs.Rdata'))

## load(file = file.path('Results','cs.Rdata'))
## class(cs$gs) <- 'logical'
## setkey(cs, gs, ord.x)

## cs[, lab:=paste0(c('gs','')[as.numeric(gs) + 1], 'RCT-',c('none','TU')[as.numeric(ord.x=='TU') + 1])]
## cs[, lab:=factor(lab)]

## pdf(file.path(figdir, 'test.pdf')
## plot(0,0, type = 'n', bty = 'n', ylim = c(0,.5), xlim = cs[.(T, 'TU'), range(ctF,ctCF)], las = 1, xlab='cases', ylab = 'power')
## cs[cf=='NTpop' & cat=='allFinalEV', points(ctF, power, pch = 15, cex = 2, col = as.numeric(lab))]
## legend('bottom', leg = cs[, unique(lab)], col = cs[, as.numeric(unique(lab))], bty = 'n', pch = 15)
## abline(v = cs[cf=='NTpop', ctCF], lty=2)
## abline(v = cs[cf=='VRpop', ctCF], lty=2)
## dev.off()

## pdf(file.path(figdir, 'test.pdf')
## plot(0,0, type = 'n', bty = 'n', ylim = c(0,.5), xlim = cs[.(T, 'TU'), range(ctF,ctCF)], las = 1, xlab='cases', ylab = 'power')
## cs[cf=='NTpop' & cat=='allFinalEV', points(ctF, power, pch = 15, cex = 2, col = as.numeric(lab))]
## legend('bottom', leg = cs[, unique(lab)], col = cs[, as.numeric(unique(lab))], bty = 'n', pch = 15)
## abline(v = cs[cf=='NTpop', ctCF], lty=2)
## abline(v = cs[cf=='VRpop', ctCF], lty=2)
## dev.off()

## break down what's happening before trial & after trial (due to roll out)
## parse out false pos/negative from correct identifications
## is there a simple linear relationship between equipoise & power? (medium risk vs high risk group, do you lose more info than equipoise gained)
## post-vacc rollout period may modulate linearity of this relationship
