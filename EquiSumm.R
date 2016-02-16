if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('login1', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/EbolaVaccSim/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T)); gc()
require(optiRum); require(RColorBrewer); require(boot); require(data.table); require(ggplot2); require(grid); require(reshape2); require(parallel)
## Simulate SWCT vs RCT vs CRCT for SL
sapply(c('multiplot.R','extractFXN.R','ggplotTheme.R'), source)

thing <- 'Equip-irsk'
## Load VaccProp & hazT
load('data/vaccProp1.Rdata')
vaccProp <- vaccProp1
vaccProp[, simNum:=1:length(vaccEff)]
## Make Figure Folder
figdir <- file.path('Figures', thing)
dir.create(figdir)

## batchdirnm <- file.path('BigResults',thing)
## fls <- list.files(batchdirnm, pattern=thing, full.names = T)
## ## fls <- fls[grepl(2305,fls)]
## eos <- extractOneSim(fls[1], indivLev=T, verbose = 0)

load(file.path('BigResults', paste0(thing, 'parmsMat','.Rdata')))

## nbtd <- parmsMat[tid==37] ##propInTrial==.05 & avHaz=='xClus' & trialStartDate=='2014-10-01'

unique(parmsMat[propInTrial==c(.05) & trialStartDate==c('2014-10-01'), list(avHaz, tid)])
nbtd <- parmsMat[tid==21,rcmdbatch]

resList <- extractSims(thing, verb=0, maxbatches=NA, nbatchDo=nbtd, indivLev = T, mc.cores=48)
resList <- procResList(resList, verb=0)

resList <- makeInfPow(resList, verb=2)
names(resList)

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
## Averted vs Spent
infAvertPow[,hazLab:=as.numeric(ihaz0Cat)]
tmp0 <- infAvertPow[avHaz=='' & propInTrial==.05 & lab %in%  c('SWCT','RCT-gs-rp')]
tmp <- rbind(tmp0[, list(trialStartDate, hazLab, lab, avHaz, type = 'avert_____', val = infAvertPC_EV)],
             tmp0[, list(trialStartDate, hazLab, lab, avHaz, type = 'spent_____', val = infSpentPC_EV)],
             tmp0[, list(trialStartDate, hazLab, lab, avHaz, type = 'avert_noEV', val = infAvertPC_noEV)],
             tmp0[, list(trialStartDate, hazLab, lab, avHaz, type = 'spent_noEV', val = infSpentPC_noEV)])
pdf(file.path(figdir, 'per capita risk averted by strata.pdf'), w = 8, h = 5)
p <- ggplot(tmp[!grepl('noEV',type)]) + geom_bar(aes(hazLab, val, fill = type), stat='identity', position='dodge') + facet_wrap(lab~trialStartDate, nc=4)  +
    guides(fill = guide_legend(keywidth = 3, keyheight = 1)) +
    ylab('per capita risk') + xlab('hazard level') + scale_color_manual(values=c('blue','red'))
print(p)
p <- ggplot(tmp[grepl('noEV',type)]) + geom_bar(aes(hazLab, val, fill = type), stat='identity', position='dodge') + facet_wrap(lab~trialStartDate, nc=4)  +
    guides(fill = guide_legend(keywidth = 3, keyheight = 1)) +
    ylab('per capita risk') + xlab('hazard level') + scale_color_manual(values=c('red','blue'))
print(p)
graphics.off()

pdf(file.path(figdir, 'per capita risk averted by strata.pdf'), w = 8, h = 5)
p <- ggplot(tmp[!grepl('noEV',type)]) + geom_bar(aes(hazLab, val, fill = type), stat='identity', position='dodge') + facet_wrap(lab~trialStartDate, nc=4)  +
    guides(fill = guide_legend(keywidth = 3, keyheight = 1)) +
    ylab('per capita risk') + xlab('hazard level') + scale_color_manual(values=c('blue','red'))
print(p)
graph
####################################################################################################
## frac of power vs infections spent 1
jpeg(file.path(figdir, paste0('infSpentPC_EV versus PoP v1.jpeg')), w = 10, h = 8, units = 'in', res = 200)
tmp <- infAvertPow[avHaz=='' & propInTrial==.1 & !ihaz0Cat %in% levels(ihaz0Cat)[c(1:4,13)] & lab!='VR']
p <- ggplot() +
    geom_point(data = tmp, aes(infSpentPC_EV, powfrac_EV, shape = lab, colour = ihaz0Cat, size = 1.5)) + 
        facet_grid(~trialStartDate) + scale_colour_manual(values=rev(brewer.pal(n=length(unique(tmp$ihaz0Cat)), "Spectral"))) +
                xlab('per capita risk spent') + ylab('fraction of information from hazard class') 
print(p)
graphics.off()
## 2
jpeg(file.path(figdir, paste0('infSpentPC_EV versus PoP v2.jpeg')), w = 12, h = 10, units = 'in', res = 200)
tmp <- infAvertPow[avHaz=='' & propInTrial==.1 & !ihaz0Cat %in% levels(ihaz0Cat)[c(1:4,13)] & lab!='VR']
p <- ggplot() +
    geom_point(data = tmp, aes(infSpentPC_EV, powfrac_EV, shape = lab, colour = lab, size = 1.5)) + 
        facet_grid(ihaz0Cat~trialStartDate) + #scale_colour_manual('Pastel1') +
                xlab('per capita risk spent') + ylab('fraction of information from hazard class') 
print(p)
graphics.off()

## 3
jpeg(file.path(figdir, paste0('infSpentPC_EV versus PoPPC v3.jpeg')), w = 10, h = 8, units = 'in', res = 200)
tmp <- infAvertPow[avHaz=='' & propInTrial==.1 & !ihaz0Cat %in% levels(ihaz0Cat)[c(1:4,13)] & lab!='VR']
p <- ggplot() +
    geom_point(data = tmp, aes(infSpentPC_EV, pow_per_infSpent_EV, shape = lab, colour = ihaz0Cat, size = 1.5)) + 
        facet_grid(~trialStartDate) + scale_colour_manual(values=rev(brewer.pal(n=length(unique(tmp$ihaz0Cat)), "Spectral"))) +
                xlab('per capita risk spent') + ylab('fraction of information from hazard class') 
print(p)
graphics.off()
####################################################################################################

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

## equipoise perturbed plot
## histogram within strata groups
## look at risk by person/strata by treatment assignment for trials
## vaccinating a greater % of people increases risk averted/spent, but still has problem of withholding treatment from individuals

##########

attach(resList)

punq <- unique(parms[,-1,with=F])
punq <- data.table(pid = 1:nrow(punq), punq)
punq[, lab:=trial]
punq[trial=='RCT' & gs==T, lab:=paste0(lab,'-gs')]
punq[trial=='RCT' & ord=='TU', lab:=paste0(lab,'-rp')]
punq[,lab:=as.factor(lab)]
setkey(punq, pid)

parms <- merge(punq, parms, by = names(punq)[!names(punq) %in% c('lab','pid')])
setkey(parms, pid, nbatch)
setcolorder(parms, c('pid','nbatch', names(parms)[!names(parms) %in% c('pid','nbatch')]))
parms
Spop <- merge(parms[,list(pid, lab, nbatch)], Spop, by = 'nbatch')
setcolorder(Spop, c('pid','lab','nbatch', names(Spop)[!names(Spop) %in% c('pid','lab', 'nbatch')]))
Spop[,assignment:=c('vacc','cont')[as.numeric(vaccDay==Inf)+1]]
setkey(Spop, Oi, pid, nbatch, simNum)
Spop[,quantile(indivRR, c(.025,.975))]

## table by individual the average infection risk in each design
irsk <- Spop[, list(.N, inf = mean(infectDay<Inf), inf_EV = mean(infectDay_EV<Inf), indivRR=unique(indivRR), Oc=Oc[1], indiv=indiv[1]), 
             list(pid,Oi,cont)]
setkey(irsk, Oi, pid)
## Average spent over all scenarios (but noting that for within-arm randomization, same individuals are always controls)
irsk[, spent := inf - inf[pid==punq[trial=='VR',pid]], list(Oi)] 
irsk[, spent_EV := inf_EV - inf[pid==punq[trial=='VR',pid]], list(Oi)]
## averted over all scenarios
irsk[, avert := inf[pid==punq[trial=='NT',pid]] - inf, list(Oi)]
irsk[, avert_EV := inf[pid==punq[trial=='NT',pid]] - inf_EV, list(Oi)]
irsk <- merge(punq[,list(pid,gs, lab)], irsk, by = 'pid')

## Average spent over worst possible randomized vaccination order assignment for non-controls
irskMax <- Spop[cont=='vacc' & vaccDay==max(vaccDay[cont=='vacc']),
                list(.N, infMax = mean(infectDay<Inf), infMax_EV = mean(infectDay_EV<Inf)),
                list(pid,Oi,vaccDay)]

setkey(irskMax, Oi, pid)
irsk <- merge(irsk, irskMax[,list(Oi,pid, infMax, infMax_EV)], by = c('Oi','pid'), all.x=T)
## Spent max
irsk[, spentMax := infMax - inf[pid==punq[trial=='VR',pid]], list(Oi)] 
irsk[, spentMax_EV := infMax_EV - inf[pid==punq[trial=='VR',pid]], list(Oi)]

##########
## Examine results
irsk[, list(inf = mean(inf), inf_EV = mean(inf_EV), infMax = mean(infMax), infMax_EV=mean(infMax_EV)), list(lab, cont)][order(lab)] ## mean infection rate without vaccine
sptmp <- irsk[, list(spent = mean(spent), spent_EV = mean(spent_EV), spentMax = mean(spentMax), spentMax_EV=mean(spentMax_EV)), list(lab, cont)][order(lab)] #
for(ii in 3:6) sptmp[[ii]] <- formatC(100*signif(sptmp[[ii]],2)) ## % risk spent on average
sptmp

cols <- punq$lab
names(cols) <- cols
cols <- c("NT"='red', "SWCT"='orange', "VR"='dark green', "RCT-gs"='dodger blue', "RCT"='dodger blue', "RCT-rp"='purple', "RCT-gs-rp" = 'purple')
## cols <- c("NT"='red', "SWCT"='orange', "VR"='dark green', "RCT-gs"='dodger blue', "RCT"='blue', "RCT-rp"='purple', "RCT-gs-rp" = 'magenta')

## histogram of risk
pdf(file.path(figdir, paste0('irsk.pdf')), w = 6.5, h = 4) #, units = 'in', res = 200)
                                        #ggplot(irsk, aes(x=inf, col = lab)) + geom_line(stat='density') + xlab('cumulative risk of infection') + 
                                        #scale_x_log10(breaks=10^c(-4:0))
                                        #    scale_x_continuous(breaks=10^c(-4:0))
ggplot(irsk[lab=='NT'], aes(x=inf)) + geom_histogram() + xlab('cumulative risk of infection') 
graphics.off()

lbrks <- c(.0001,.0005,.001,.005,.01,.05,.1,.5,1)

## histogram of risk spent
pdf(file.path(figdir, paste0('irsk spent hist.pdf')), w = 6.5, h = 4)#, units = 'in', res = 200)
ggplot(irsk[!lab %in% c('VR','NT')], aes(x=spent_EV, fill=lab)) + geom_histogram() + facet_wrap(~lab, ncol=1) + xlim(-.04, .08)
graphics.off()

## density lines risk spent
pdf(file.path(figdir, paste0('irsk spent.pdf')), w = 6.5, h = 4)#, units = 'in', res = 200)
adj <- 2
p <- ggplot() + 
    geom_line(data = irsk[!lab %in% c('NT','VR')], aes(x=spent_EV, col=lab, linetype = gs), stat='density', adjust=adj) +
              scale_color_manual(values=cols) + xlim(.01, .3) + ylim(0,32) + xlab('per capita infection risk spent')  + ggtitle('average risk spent per individual')
print(p) ## print(p+scale_x_log10(breaks=lbrks))
## spentMax_EV for controls in non-rp RCTs, spent_EV for controls in rp RCTs, spentMax_EV for anyone in SWCT
p <- ggplot() + 
    geom_line(data = irsk[lab %in% c('RCT','RCT-gs')], aes(x=spentMax_EV, col=lab, linetype = gs), stat='density', adjust=adj) +
    geom_line(data = irsk[lab %in% c('RCT-rp','RCT-gs-rp') & cont=='cont'], aes(x=spent_EV, col=lab, linetype = gs), stat='density', adjust=adj) + 
    geom_line(data = irsk[lab %in% c('SWCT')], aes(x=spentMax_EV, col=lab), stat='density', adjust=adj) +                       
              scale_color_manual(values=cols) + xlim(.01, .3) + ylim(0,32) + xlab('per capita infection risk spent')  + ggtitle('maximum risk spent per individual')
print(p) ## print(p+scale_x_log10(breaks=lbrks))
## spentMax_EV
graphics.off()

irsk[,lab:=factor(lab, levels(lab)[c(1,3,2,4,5,6)])]
irsk <- irsk[order(Oc,indivRR,lab)]

## spent
jpeg(file.path(figdir, paste0('irsk spent bars.jpeg')), w = 6.5, h = 4, units = 'in', res = 200)
ggplot(irsk[!lab %in% c('VR','NT')], aes(x=indiv, y=spent, fill=cont)) + 
    geom_bar(stat='identity') + facet_wrap(~lab) + ylim(-.05,.2)
graphics.off()

jpeg(file.path(figdir, paste0('irsk spentEV bars.jpeg')), w = 6.5, h = 4, units = 'in', res = 200)
ggplot(irsk[!lab %in% c('VR','NT')], aes(x=indiv, y=spent_EV, fill=cont)) + 
    geom_bar(stat='identity') + facet_wrap(~lab) + ylim(-.05,.2)
graphics.off()

## avert
jpeg(file.path(figdir, paste0('irsk avert bars.jpeg')), w = 6.5, h = 4, units = 'in', res = 200)
ggplot(irsk[!lab %in% c('VR','NT')], aes(x=indiv, y=avert, fill=cont)) + 
    geom_bar(stat='identity') + facet_wrap(~lab) + ylim(-.05,.2)
graphics.off()

jpeg(file.path(figdir, paste0('irsk avertEV bars.jpeg')), w = 6.5, h = 4, units = 'in', res = 200)
ggplot(irsk[!lab %in% c('VR','NT')], aes(x=indiv, y=avert_EV, fill=cont)) + 
    geom_bar(stat='identity') + facet_wrap(~lab) + ylim(-.05,.2)
graphics.off()

## inf
jpeg(file.path(figdir, paste0('irsk inf bars.jpeg')), w = 6.5, h = 4, units = 'in', res = 200)
ggplot(irsk[lab=='NT'], aes(x=indiv, y=inf, fill=cont)) + 
    geom_bar(stat='identity') + ylim(0,1)
graphics.off()


