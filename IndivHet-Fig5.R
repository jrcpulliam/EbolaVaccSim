
## infection risk
irsk[,cluster:=factor(cluster)]
p <- ggplot(irsk[lab=='NT'], aes(x=ordShow, y=inf, fill=cluster)) +  ylab('cumulative infection risk') + xlab('individual') +
    geom_bar(stat='identity', width=1) + theme(legend.key.size = unit(.1, "cm")) + ggtitle('infection risk without vaccination') +
        theme(legend.position="right") #+ theme(axis.title.y = element_text(angle=0))
ggsave(file.path(figdir, paste0('irsk inf bars.pdf')), plot=p, w=wid, h=heig, units='in')


## Conditional on arms & order randomization
tmp <- irsk[(arm==armShown & type=='cond' &  grepl('RCT',lab)) | (type=='condvd' & lab=='SWCT' & exmpl==T)]
tmp <- tmp[lab!='RCT-rp']
####################################################################################################
## spent
ylim <- tmp[,range(spent,spent_EV)]
p <- ggplot(tmp, aes(x=ordShowArm, y=spent, fill=armShown)) + ggtitle('risk spent, conditional on arm without EV') +
    geom_bar(stat='identity', width=1) + facet_grid(lab ~ .) +  ylab('risk') + xlab('individual') + ylim(ylim[1],ylim[2])
ggsave(file.path(figdir, paste0('irsk spent bars.jpeg')), p, w = wid, h = heig, units = 'in')


p <- ggplot(tmp, aes(x=ordShowArm, y=spent_EV, fill=armShown)) + ggtitle('risk spent, conditional on arm with EV') +
    geom_bar(stat='identity', width=1) + facet_grid(lab ~ .) +  ylab('risk')+ xlab('individual') + ylim(ylim[1],ylim[2])
ggsave(file.path(figdir, paste0('irsk spentEV bars.jpeg')), p, w = wid, h = heig, units = 'in')

## avert
ylim <- tmp[,range(avert,avert_EV)]
p <- ggplot(tmp, aes(x=ordShowArm, y=avert, fill=armShown)) + ggtitle('risk averted, conditional on arm without EV') +
    geom_bar(stat='identity', width=1) + facet_grid(lab ~ .) +  ylab('risk')+ xlab('individual') + ylim(ylim[1],ylim[2])
ggsave(file.path(figdir, paste0('irsk avert bars.jpeg')), p, w = wid, h = heig, units = 'in')


p <- ggplot(tmp, aes(x=ordShowArm, y=avert_EV, fill=armShown)) + ggtitle('risk averted, conditional on arm with EV') +
    geom_bar(stat='identity', width=1) + facet_grid(lab ~ .) +  ylab('risk')+ xlab('individual') + ylim(ylim[1],ylim[2])
ggsave(file.path(figdir, paste0('irsk avertEV bars.jpeg')), p, w = wid, h = heig, units = 'in')
 
## both
ylim <- tmp[,range(-avert_EV,spent_EV, -avert, spent)]

p <- ggplot(tmp) + ggtitle('risk spent, conditional on arm without EV') +
    geom_hline(yintercept=.05) +
        geom_bar(aes(x=ordShowArm, y=spent, fill=armShown), stat='identity', width=1) +
            geom_bar(aes(x=ordShowArm, y=-avert, fill=armShown), stat='identity', width=1, alpha = .8) +
                facet_grid(lab ~ .) +  ylab('risk') + xlab('individual') + ylim(ylim[1],ylim[2])
ggsave(file.path(figdir, paste0('irsk spent & avert bars.jpeg')), p, w = wid, h = heig, units = 'in')


p <- ggplot(tmp) + ggtitle('risk spent, conditional on arm with EV') +
    geom_hline(yintercept=.05) +
        geom_bar(aes(x=ordShowArm, y=spent_EV, fill=armShown), stat='identity', width=1) +
            geom_bar(aes(x=ordShowArm, y=-avert_EV, fill=armShown), stat='identity', width=1, alpha = .8) +
                facet_grid(lab ~ .) +  ylab('risk') + xlab('individual') + ylim(ylim[1],ylim[2])
ggsave(file.path(figdir, paste0('irsk spent & avert EV bars.jpeg')), p, w = wid, h = heig, units = 'in')

## SB version
ylim <- tmp[,range(-avert_EV,spent_EV, -avert, spent)]
tmp[, cols:=armShown]; tmp[cols=='cont',cols:='red']; tmp[cols=='vacc',cols:='dodger blue']
pdf(file.path(figdir, paste0('irsk spent & avert SB.pdf')), w = wid, h = heig)
par(mfrow=c(4,1), mar = c(0,3,1,0), oma = c(1,1,0,0))
for(ll in tmp[,unique(lab)]) {
    with(tmp[lab==ll], plot(ordShowArm, spent_EV, xlab='individual', ylab='risk spent', bty = 'n', type = 'h', col = cols, ylim = ylim, las = 1, xaxt='n', main =ll))
    with(tmp[lab==ll], points(ordShowArm, -avert_EV, type = 'h', col = makeTransparent(cols, alpha = 250)))
    abline(h=0, lty = 1)
    abline(h=.05, lty = 2, col='gray')    
}
title(xlab='individual',outer=T)
title(ylab='risk',outer=T)
graphics.off()


## SB version short
ylim <- tmp[,range(-avert_EV,spent_EV, -avert, spent)]
tmp[, cols:=armShown]; tmp[cols=='cont',cols:='red']; tmp[cols=='vacc',cols:='dodger blue']
pdf(file.path(figdir, paste0('irsk spent & avert SB.pdf')), w = wid, h = heig)
par(mfrow=c(2,1), mar = c(0,5,1,0), oma = c(1,2,0,0))
ll='RCT-gs-rp'
with(tmp[as.numeric(cluster) <= 8 & lab==ll], plot(ordShowArm, spent, xlab='individual', ylab='spent', bty = 'n', type = 'h', col = cols, ylim = c(0,.2), las = 1, xaxt='n', main =''))
abline(h=0, lty = 1)
abline(h=.05, lty = 2, col='dark gray', lwd  =2)
legend('topright', leg = c('cont','vacc'), col = c('red','dodger blue'), pch = 15, bty = 'n')
with(tmp[as.numeric(cluster) <= 8 & lab==ll], plot(ordShowArm, -avert, xlab='individual', ylab='averted', bty = 'n', type = 'h', col = cols, ylim = c(-.2,0), las = 1, xaxt='n', main =''))
title(ylab='risk',outer=T, line = 0)
graphics.off()


## SB risk
tmp <- irsk[lab=='NT']
ylim <- tmp[,range(inf)]
pdf(file.path(figdir, paste0('irsk inf risk SB.pdf')), w = wid, h = heig)
par(mfrow=c(2,1), mar = c(0,5,1,0), oma = c(1,2,0,0))
with(tmp[as.numeric(cluster) <= 8], plot(ordShow, inf, xlab='individual', ylab='', bty = 'n', type = 'h', col = rainbow(8)[as.numeric(cluster)], ylim = c(0,.2), las = 1, xaxt='n', main =''))
title(ylab='risk',outer=T, line = 0)
graphics.off()

####################################################################################################

## marginal on arms & order randomization
cols <- c('w/o EV' = 'dodger blue', 'w/ EV' = 'purple')
tmp <- irsk[type=='marg' & !lab %in% c('NT','VR')]
####################################################################################################
## spent
ylim <- tmp[,range(spent,spent_EV)]

p <- ggplot(tmp) + ggtitle('risk spent, marginal on randomization') +
        geom_bar(aes(x=ordShow, y=spent, fill=ev), data.table(tmp, ev = 'w/o EV'), stat='identity', width=1) +
        geom_bar(aes(x=ordShow, y=spent_EV, fill=ev), data.table(tmp[lab!='SWCT'], ev = 'w/ EV'), stat='identity', width=1) +
        scale_fill_manual(values=cols) +
        facet_grid(lab ~ .) +  ylab('risk') + xlab('individual') + ylim(ylim[1],ylim[2])
ggsave(file.path(figdir, paste0('irsk spent bars MARG.jpeg')), p, w = wid, h = heig, units = 'in')

## avert
ylim <- c(0, tmp[,max(avert,avert_EV)])

p <- ggplot(tmp) + ggtitle('risk averted, marginal on randomization') +
    geom_bar(aes(x=ordShow, y=avert_EV, fill=ev), data.table(tmp[lab!='SWCT'], ev = 'w/ EV'), stat='identity', width=1) +
    geom_bar(aes(x=ordShow, y=avert, fill=ev), data.table(tmp, ev = 'w/o EV'), stat='identity', width=1) +
        scale_fill_manual(values=cols) +
            facet_grid(lab ~ .) +
                ylab('risk')+ xlab('individual') + ylim(ylim[1],ylim[2])
ggsave(file.path(figdir, paste0('irsk avert bars MARG.jpeg')), p, w = wid, h = heig, units = 'in')


p <- ggplot(SpopH[simNum==1 & nbatch==3073], aes(Date, clusHaz*10^5, group = OcOrd, col=OcOrd)) + geom_line() + theme(legend.position="top") +
     theme(legend.key.size = unit(.1, "cm")) + ylab('daily infection hazard (per 100,000)') + ggtitle('mean cluster hazard trends')
ggsave(file.path(figdir, paste0('haz traj.jpeg')), p, w = wid, h = heig, units = 'in')

## **check why # of SWCT with posv is not the same as other designs** could be just running different sims?


p <- ggplot(triSumm, aes(above_EV, power, col = lab)) + geom_point() + ylim(0,1) + xlab(paste0('expected # subjects spending >', threshold, ' infection risk'))
ggsave(file.path(figdir, paste0('pow trhes.jpeg')), p, w = wid, h = heig, units = 'in')


## add error bars to inf spent & frac of information
## try fraction of information per person on y-axis
## show inf spent/averted on same plot

## equipoise perturbed plot
## histogram within strata groups
## look at risk by person/strata by treatment assignment for trials
## vaccinating a greater % of people increases risk averted/spent, but still has problem of withholding treatment from individuals
