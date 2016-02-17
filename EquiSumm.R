if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('login1', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/EbolaVaccSim/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T)); gc()
require(optiRum); require(RColorBrewer); require(boot); require(data.table); require(ggplot2); require(grid); require(reshape2); require(parallel)
moveFront <- function(dt, fnames) setcolorder(dt, c(fnames, names(dt)[!names(dt) %in% fnames]))
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

## resList <- makeInfPow(resList, verb=2)
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
punq[,lab:=factor(lab, levels =levels(lab)[c(1:3,5,4,6:7)])]
setkey(punq, pid)

parms <- merge(punq, parms, by = names(punq)[!names(punq) %in% c('lab','pid')])
setkey(parms, pid, nbatch)
setcolorder(parms, c('pid','nbatch', names(parms)[!names(parms) %in% c('pid','nbatch')]))
parms
Spop <- merge(parms[,list(pid, lab, nbatch)], Spop, by = 'nbatch')
setcolorder(Spop, c('pid','lab','nbatch', names(Spop)[!names(Spop) %in% c('pid','lab', 'nbatch')]))
Spop[,arm:=c('vacc','cont')[as.numeric(vaccDay==Inf)+1]]
setkey(Spop, Oi, pid, nbatch, simNum)
Spop[,quantile(indivRR, c(.025,.975))]

## table by individual the average infection risk in each design
irskMarg <- Spop[, list(.N, inf = mean(infectDay<Inf), inf_EV = mean(infectDay_EV<Inf) , indivRR=unique(indivRR), Oc=Oc[1], type = 'marg', arm = NA), 
             list(pid,Oi)]
irskCond <- Spop[!lab %in% c('VR','NT','SWCT'), ## conditional on control/vacc randomization assignment
             list(.N, inf = mean(infectDay<Inf), inf_EV = mean(infectDay_EV<Inf), indivRR=unique(indivRR), Oc=Oc[1], type = 'cond'), 
             list(pid,Oi,arm)]
## maximum spent given worst possible vaccination order (only for random ordered trials)
irskMax <- Spop[!lab %in% c('VR','NT') & arm=='vacc' & !grepl('rp',lab) & vaccDay==max(vaccDay[arm=='vacc']),
                list(.N, inf = mean(infectDay<Inf), inf_EV = mean(infectDay_EV<Inf), indivRR=unique(indivRR), Oc=Oc[1], arm='vacc', type='max'),
                list(pid,Oi,vaccDay)]
## conditional on exact vaccination day (note borrowing control info across all vaccDay==Inf, i.e. whether the cluster they're in is vaccinated early/late)
irskCondvd <- Spop[!lab %in% c('VR','NT'),
                list(.N, inf = mean(infectDay<Inf), inf_EV = mean(infectDay_EV<Inf), indivRR=unique(indivRR), Oc=Oc[1], type='condvd'),
                list(pid,Oi,arm,vaccDay)]

irsk <- rbindlist(list(irskMarg, irskCond, irskMax, irskCondvd), use.names=T, fill=T)
setkey(irsk, Oi, pid)
irsk <- merge(punq[,list(pid,gs, lab)], irsk, by = 'pid')
moveFront(irsk, c('pid','lab','gs','type','arm','vaccDay'))
irsk
irsk[, unique(type), lab]
irsk[Oi==1]

## Average risk spent versus marginal VR 
irsk[, spent   := inf    - inf[lab=='VR'], list(Oi)]
irsk[, spent_EV:= inf_EV - inf[lab=='VR'], list(Oi)]
## Averted risk versus marginal NT
irsk[, avert   :=    inf[lab=='NT'] - inf, list(Oi)]
irsk[, avert_EV:= inf[lab=='NT'] - inf_EV, list(Oi)]
setkey(irsk, Oi, pid)

irsk[Oi==1]

##########
## Examine results
irsk[type!='condvd', list(inf = mean(inf), inf_EV = mean(inf_EV)), list(lab, arm, type)][order(lab)] ## mean infection rate without vaccine
sptmp <- irsk[type!='condvd', list(spent = mean(spent), spent_EV = mean(spent_EV)), list(lab, arm, type)][order(lab)] #
for(ii in 4:5) sptmp[[ii]] <- formatC(100*signif(sptmp[[ii]],2)) ## % risk spent on average
sptmp

lbrks <- c(.0001,.0005,.001,.005,.01,.05,.1,.5,1)
cols <- punq$lab
names(cols) <- cols
cols <- c("NT"='red', "SWCT"='orange', "VR"='dark green', "RCT-gs"='dodger blue', "RCT"='dodger blue', "RCT-rp"='purple', "RCT-gs-rp" = 'purple')
## cols <- c("NT"='red', "SWCT"='orange', "VR"='dark green', "RCT-gs"='dodger blue', "RCT"='blue', "RCT-rp"='purple', "RCT-gs-rp" = 'magenta')

## get clusters in risk-order
cOrd <- irsk[lab=='NT', list(inf=mean(inf)), Oc][order(-inf)]
cOrd[,cluster:=1:nrow(cOrd)]
irsk$cluster <- NULL
irsk <- merge(irsk, cOrd[,list(Oc, cluster)], by = 'Oc')
irsk[,cluster:=factor(cluster)]

## order individuals within clusters by risk for ease of display
iord <- irsk[lab=='NT',list(Oi,cluster,inf)][order(cluster,-inf)]
iord[,ordShow:=1:6000]
irsk$ordShow <- NULL
irsk <- merge(irsk, iord[,list(Oi,ordShow)], by = 'Oi')

## what arm to label individuals as for RCTs?
irsk[, armShown:=arm]
irsk[grepl('RCT',lab) & as.numeric((Oi-1) %% (300) < 150), armShown:='cont']
irsk[grepl('RCT',lab) & as.numeric((Oi-1) %% (300) >= 150), armShown:='vacc']

irsk0 <- copy(irsk)
irsk <- copy(irsk0)

## order individuals within clusters & armShown by risk for ease of display
iord <- irsk[lab=='NT',list(Oi,cluster,inf)]
iord <- merge(iord, unique(irsk[lab=='RCT' & arm==armShown, list(Oi, armShown)]))
iord <- iord[order(cluster,armShown, -inf)]
iord[,ordShowArm:=1:6000]
irsk$ordShowArm <- NULL
irsk <- merge(irsk, iord[,list(Oi,ordShowArm)], by = 'Oi')

irsk <- irsk[order(Oc,indivRR,lab)]
irsk[type=='cond' &  !lab %in% c('VR','NT') & pid==2 & Oi==1]

## need exemplar SWCT: pick arbitrary example of cluster-day assignment to display
## **(later do probably do exemplars for all of them with exact same randomization throughout)
clusVD <- irsk[(arm==armShown & type=='condvd' &  grepl('SWCT',lab)),list(Oc = 1:20, vaccDay=unique(vaccDay))]
setkey(clusVD, Oc, vaccDay)
irsk[,exmpl:=F]

setkey(irsk, Oc, vaccDay)
irsk[clusVD][type=='condvd' & 'SWCT'==lab][["exmpl"]] <- rep(T,6000)
irsk[clusVD][type=='condvd' & 'SWCT'==lab]
setkey(irsk, Oi, pid)

## histogram of risk
pdf(file.path(figdir, paste0('irsk.pdf')), w = 6.5, h = 4) 
ggplot(irsk[lab=='NT'], aes(x=inf)) + geom_histogram() + xlab('cumulative risk of infection') 
graphics.off()

## density lines risk spent
pdf(file.path(figdir, paste0('irsk spent.pdf')), w = 6.5, h = 4)#, units = 'in', res = 200)
adj <- 1.5
p <- ggplot() + 
    geom_line(data = irsk[type=='marg'], aes(x=spent_EV, col=lab, linetype = gs), stat='density', adjust=adj) +
              scale_color_manual(values=cols) + xlim(.005, .3) + ylim(0,32) + xlab('per capita infection risk spent')  + ggtitle('average risk spent per individual')
print(p) ## print(p+scale_x_log10(breaks=lbrks))
p <- ggplot() + 
    geom_line(data = irsk[(type=='max' & lab=='SWCT') | (type='cond' & arm=='cont' & lab=='RCT-gs-rp')],
              aes(x=spent_EV, col=lab, linetype = gs), stat='density', adjust=adj) +
              scale_color_manual(values=cols) + xlim(.005, .3) + ylim(0,32) + xlab('per capita infection risk spent')  + ggtitle('max risk spent per individual')
print(p) ## print(p+scale_x_log10(breaks=lbrks))
graphics.off()

tmp <- irsk[(arm==armShown & type=='cond' &  grepl('RCT',lab)) | (type=='condvd' & lab=='SWCT' & exmpl==T)]
res <- 300

####################################################################################################
## spent
jpeg(file.path(figdir, paste0('irsk spent bars.jpeg')), w = 6.5, h = 4, units = 'in', res = res)
ggplot(tmp, aes(x=ordShowArm, y=spent, fill=armShown)) + ggtitle('risk spent, conditional on arm without EV') +
    geom_bar(stat='identity') + facet_wrap(~lab, ncol=2) + ylim(-.05,.2) + ylab('risk') + xlab('individual')
graphics.off()

jpeg(file.path(figdir, paste0('irsk spentEV bars.jpeg')), w = 6.5, h = 4, units = 'in', res = res)
ggplot(tmp, aes(x=ordShowArm, y=spent_EV, fill=armShown)) + ggtitle('risk spent, conditional on arm with EV') +
    geom_bar(stat='identity') + facet_wrap(~lab, ncol=2) + ylim(-.05,.2) + ylab('risk')+ xlab('individual')
graphics.off()

## avert
jpeg(file.path(figdir, paste0('irsk avert bars.jpeg')), w = 6.5, h = 4, units = 'in', res = res)
ggplot(tmp, aes(x=ordShowArm, y=avert, fill=armShown)) + ggtitle('risk averted, conditional on arm without EV') +
    geom_bar(stat='identity') + facet_wrap(~lab, ncol=2) + ylim(-.05,.2) + ylab('risk')+ xlab('individual')
graphics.off()

jpeg(file.path(figdir, paste0('irsk avertEV bars.jpeg')), w = 6.5, h = 4, units = 'in', res = res)
ggplot(tmp, aes(x=ordShowArm, y=avert_EV, fill=armShown)) + ggtitle('risk averted, conditional on arm with EV') +
    geom_bar(stat='identity') + facet_wrap(~lab, ncol=2) + ylim(-.05,.2) + ylab('risk')+ xlab('individual')
graphics.off()
####################################################################################################

## inf

jpeg(file.path(figdir, paste0('irsk inf bars.jpeg')), w = 6.5, h = 4, units = 'in', res = res)
ggplot(irsk[lab=='NT'], aes(x=ordShow, y=inf, col=factor(cluster))) +  ylab('cumulative infection risk \n(over year after trial start)') + xlab('individual') +
    geom_bar(stat='identity') + ylim(0,1) + theme(legend.key.size = unit(.1, "cm")) + ggtitle('infection risk without vaccination') +
        theme(legend.position="top")
graphics.off()

##
SpopH[,cluster:=factor(cluster)]
cOrd[,OcOrd:=factor(cluster)]

SpopH$OcOrd <- NULL
SpopH <- merge(SpopH, cOrd[,list(Oc, OcOrd)], by = 'Oc')

SpopH[simNum==1 & nbatch==3073]

pdf(file.path(figdir, paste0('haz traj.pdf')), w = 6.5, h = 4)#, units = 'in', res = res)
ggplot(SpopH[simNum==1 & nbatch==3073], aes(Date, clusHaz*10^5, group = OcOrd, col=OcOrd)) + geom_line() + theme(legend.position="top") +
     theme(legend.key.size = unit(.1, "cm")) + ylab('daily infection hazard (per 100,000)') + ggtitle('mean cluster hazard trends')
graphics.off()

