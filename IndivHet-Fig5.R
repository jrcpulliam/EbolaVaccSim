if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('sbellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('nid', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/EbolaVaccSim/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T)); gc()
require(optiRum); require(RColorBrewer); require(boot); require(data.table); require(ggplot2); require(grid); require(reshape2); require(parallel)
sapply(c('multiplot.R','extractFXN.R','ggplotTheme.R'), source)
wid <- 6.5
heig <- 4
res <- 300

thing <- 'Equip-Fig5-v4'
figdir <- file.path('Figures', thing)
dir.create(figdir)
fls <- list.files('BigResults', pattern = paste0(thing,'-'), full.names=T)
colsv <- c("NT"='red', "SWCT"='orange', "VR"='dark green', "RCT-gs"='dodger blue', "RCT"='light blue', "RCT-rp"='purple', "RCT-gs-rp" = 'purple')

load(file.path('BigResults', paste0(thing, 'parmsMat','.Rdata')))
tvars <- c("trialStartDate", "propInTrial", "avHaz", "indivRRSeed" , "HazTrajSeed",'numClus','clusSuze','hazType','nbsize','mu','cvClus','cvClusTime','sdLogIndiv','weeklyDecay','cvWeeklyDecay','hazIntUnit')
tvars <- tvars[tvars %in% colnames(parmsMat)]
tpop <- unique(parmsMat[,tvars, with=F]) ## make sure to add other variables unique to population
tpop$tid <- 1:nrow(tpop)

irsk <- punq <- data.table()
for(ii in 1:length(fls)) {
    load(fls[ii])
    tmp <- merge(resList$punq,resList$SpopWorst,by = c('pid', 'propInTrial','lab'))
    tmp <- merge(tmp, resList$irsk[type=='marg', list(caseSpent = 6000*mean(spent_EV)), pid], by = 'pid')
    itmp <- resList$irsk
    itmp$tid <- ii ## **careful here**, need to make sure it always matches flow from parmsMat
    irsk <- rbind(irsk, itmp)
    speedpow <- with(resList, merge(parms[, list(pid, nbatch)], finTrials[,list(tcal, vaccEff,vaccGood, cvr,nbatch)], by = 'nbatch'))
    tmp <- merge(tmp, speedpow[, list(tcal = mean(tcal)), pid], by = 'pid')
    punq <- rbind(punq, tmp)
}
punq[,trialStartDate:=as.Date(trialStartDate)]
punq[,date:=format.Date(trialStartDate, '%b-%y')]
plunq <- punq[threshold==.05, list(lab, power, trialStartDate, threshold, above, above_EV, caseSpent, totCase, totCase_EV,avHaz, tcal,date)]

## Create plot that shows risk averted

####################################################################################################
## infection risk
itmp <- irsk[tid==1]

p <- ggplot(itmp[lab=='NT'], aes(x=ordShow, y=inf, fill=cluster)) +  ylab('cumulative infection risk') + xlab('individual') +
    geom_bar(stat='identity', width=1) + theme(legend.key.size = unit(.1, "cm")) + ggtitle('infection risk without vaccination') +
        theme(legend.position="right") #+ theme(axis.title.y = element_text(angle=0))
print(p)
ggsave(file.path(figdir, paste0('itmp inf bars.pdf')), plot=p, w=wid, h=heig, units='in')

## Conditional on arms & order randomization
itmp <- itmp[((arm==armShown & type=='cond' &  grepl('RCT',lab)) | (type=='condvd' & lab=='SWCT' & exmpl==T))]
itmp[, cols:=armShown]; itmp[cols=='cont',cols:='red']; itmp[cols=='vacc',cols:='dodger blue']
itmp <- itmp[lab!='RCT-rp']

####################################################################################################
## SB version (not gg plot)
## Look at individuals infection risk, but dividing them into their randomization structure (for comparison to below)
xlim <- c(0,6000)
par(mfrow=c(6,1), mar = c(0,3,1,0), oma = c(1,1,0,0))
for(ll in itmp[,unique(lab)]) {
    itmp[lab==ll, plot(ordShowArm, inf, xlab='individual', ylab='risk spent', bty = 'n', type = 'h', col = cols, ylim = c(0,1), xlim=xlim, 
             las = 1, xaxt='n', main =ll)]
    itmp[lab==ll & Oi %%300==1, text(ordShowArm, .9, cluster)]
}


ylim <- itmp[,range(-avert_EV,spent_EV, -avert, spent)]
xlim <- c(0,6000)
##pdf(file.path(figdir, paste0('irsk spent & avert SB.pdf')), w = wid, h = heig)
par(mfrow=c(7,1), mar = c(0,3,1,0), oma = c(1,1,0,0))
itmp[pid==1, plot(Oi, inf, xlab='individual', ylab='risk spent', bty = 'n', type = 'h', col = cols, ylim = c(0,1), xlim=xlim, 
         las = 1, xaxt='n', main =ll)]
itmp[pid==1 & Oi %%300==1, text(Oi, .9, cluster)]
for(ll in itmp[,unique(lab)]) {
    itmp[lab==ll, plot(Oi, spent_EV, xlab='individual', ylab='risk spent', bty = 'n', type = 'h', col = cols, ylim = ylim, xlim=xlim, 
las = 1, xaxt='n', main =ll)]
    with(itmp[lab==ll], points(Oi, -avert_EV, type = 'h', col = makeTransparent(cols, alpha = 250)))
    abline(h=0, lty = 1)
    abline(h=.05, lty = 2, col='gray')    
    itmp[lab==ll & Oi %%300==1, text(Oi, ylim[2]*.9, cluster)]
}
title(xlab='individual',outer=T)
title(ylab='risk',outer=T)
##graphics.off()

itmp[cluster%in% c(1,6), list(maxinf = max(inf), maxspent=max(spent_EV), maxavert=max(avert_EV)), .(cluster,Oc)]

## so the problems are that 1) it's unlikely i'll be able to see the details for 6000, & 2) for all the trials. If I just show the average by risk group within each cluster, that shows complete info except doesn't highlight the breakdown by risk group.

## the goal is to show how each design affects risk spent by specific groups differently

names(itmp)
itmp[lab=='RCT-gs-rp-cvd' & arm=='cont', unique(cluster)]
itmp[lab=='RCT-gs-rp-cvd' & cluster==2, unique(arm)]
Spop <- resList$Spop
Spop[lab=='RCT-gs-rp-cvd' & arm=='contVD', unique(cluster)]
Spop[lab=='RCT-gs-rp-cvd' & simNum==1 & sim==1 & indiv %% 150 %in% c(1,0), .(lab, Oi, Oc, indiv, cluster, arm)][order(indiv)]
## for some reason only first cluster has a group labeled as contVD for these trials


mtmp <- irsk[type=='cond', .(inf=mean(inf), spent=mean(spent), spent_EV=mean(spent_EV), avert_EV=mean(avert_EV), avert=mean(avert), n = length(avert_EV), Oi=Oi, Oc=Oc, arm=arm), .(cluster, armO, indivRR, lab)][order(cluster, armO, -indivRR)]
mtmp[,ord:=order(Oc, indivRR, Oi)]
mtmp[, cols:=armO]; mtmp[grepl('cont',cols),cols:='red']; mtmp[grepl('vacc',cols),cols:='dodger blue']
mtmp[grepl('Excl',arm)]
mtmp[grepl('Excl',arm) & armO=='cont'] ## not seeing controls in the vaccexcluded group
## also not seeing the averted risk amongst the maxRR plot


par(mfrow=c(6,1), mar = c(0,3,1,0), oma = c(1,1,0,0))
for(ll in mtmp[,unique(lab)]) {
    mtmp[lab==ll, plot(ord, spent_EV, xlab='individual', ylab='risk spent', bty = 'n', type = 'h', col = cols, 
             las = 1, xaxt='n', main =ll)]
    with(mtmp[lab==ll], points(ord, -avert_EV, type = 'h', col = makeTransparent(cols, alpha = 250)))
    abline(h=0, lty = 1)
}
##graphics.off()

## look at all yvars for each trial type to see if things look right
pdf(file.path(figdir, paste0('irsk InfSpentAvert by trial.pdf')), w = wid, h = heig)
for(ll in mtmp[,unique(lab)]) {
    par(mfrow=c(6,1), mar = c(0,3,1,0), oma = c(1,1,0,0))
    for(yvar in c('inf', 'avert','avert_EV', 'spent', 'spent_EV')) {
        ymax <- ifelse(yvar=='inf', 1, .15)
        mtmp[lab==ll, plot(ord, get(yvar), xlab='individual', ylab='risk spent', bty = 'n', type = 'h', col = cols, 
                 ylim = c(0, ymax)
                 , las = 1, xaxt='n', main =yvar)]
    abline(h=0, lty = 1)
        mtext(ll, outer = T, side = 1, line = -1)
    }
}
graphics.off()
## OK seems like I'm missing the m


with(itmp[lab==ll], plot(ordShow, spent_EV, xlab='individual', ylab='risk spent', bty = 'n', type = 'h', col = cols, ylim = ylim, las = 1, xaxt='n', main =ll, xlim = c(0,600)))

## risk averted & avertable risk not averted for each design
## avert_EV & spent_EV

## add error bars to inf spent & frac of information
## try fraction of information per person on y-axis
## show inf spent/averted on same plot

## equipoise perturbed plot
## histogram within strata groups
## look at risk by person/strata by treatment assignment for trials
## vaccinating a greater % of people increases risk averted/spent, but still has problem of withholding treatment from individuals

par(mfrow=c(3,1))
xmax <- 600
itmp[lab=='RCT',plot(ordShowArm, inf, type = 'h', col = cols, xlim = c(0,xmax))]
itmp[lab=='RCT',plot(ordShowArm, spent_EV, type = 'h', col = cols, xlim = c(0,xmax))]
itmp[lab=='RCT',plot(ordShowArm, avert_EV, type = 'h', col = cols, xlim = c(0,xmax))]
