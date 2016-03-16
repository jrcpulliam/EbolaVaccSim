
####################################################################################################
## Analytical equipoise functions

require(data.table); require(ggplot2); require(scales)
eqd <- as.data.table(expand.grid(infRisk = c(.005, .01, .05, .1)
                               , probVaccWorks = seq(.1, .7, by = .1)
                               , efficacy = seq(.6, .9, by = .1)
                               , probSAE = 1/c(200,500,1000,10^4,10^5)
                               , cfr = c(.1,.3,.5,.7,.9)))#seq(0, 1, by = .1)))

eqd[, vaccCFR := probSAE + (1-probSAE) * (probVaccWorks * (infRisk * (1-efficacy) * cfr) + (1-probVaccWorks) * infRisk * cfr)]
eqd[, contCFR := infRisk * cfr]
eqd[, excessCFR := contCFR - vaccCFR]
eqd <- eqd[rev(order(cfr))]

eqd[infRisk==.05 & probVaccWorks==.5 & efficacy==.8 & cfr== .7]

cols <- colorRampPalette(c('blue','purple','red'))(length(unique(eqd$cfr)))

pdf('~/Desktop/eq.pdf', w =6.5, h = 4)
for(ps in c(1/200, 10^-4)) { 
    tmp <- eqd[probSAE==ps]
    p <- ggplot(tmp) + facet_grid(efficacy~infRisk) + scale_color_manual(values=cols) + # palette='RdBu') + 
        geom_line(aes(probVaccWorks, excessCFR, col = factor(cfr), group=cfr))  + 
            geom_hline(yintercept=0, linetype = 2, col = 'black') + guides(col = guide_legend(rev=T, title='CFR') ) +
                ylab('excess risk of death in control arm') + xlab('probability vaccine works')  + 
                    ggtitle(paste0('pSAE = ',ps)) + ylim(-.01, .06)
    print(p)
}
graphics.off()



pdf('~/Desktop/eq2.pdf', w =6.5, h = 4)
tmp <- eqd[efficacy==.8 & probSAE>10^-5]
p <- ggplot(tmp) + facet_grid(probSAE~infRisk) + scale_color_manual(values=cols) + # palette='RdBu') + 
    geom_line(aes(probVaccWorks, excessCFR, col = factor(cfr), group=cfr))  + 
        geom_hline(yintercept=0, linetype = 2, col = 'black') + guides(col = guide_legend(rev=T, title='CFR') ) +
            ylab('excess risk of death in control arm') + xlab('probability vaccine works')  + 
                ylim(-.01, .06) + scale_x_continuous(breaks=seq(0,1,b=.2))
print(p)
graphics.off()

pdf('~/Desktop/eq3.pdf', w =6.5, h = 4)
tmp <- eqd[efficacy==.8 & probSAE==10^-4 & infRisk %in% c(.01,.05)]
p <- ggplot(tmp) + facet_grid(~infRisk) + scale_color_manual(values=cols) + # palette='RdBu') + 
    geom_line(aes(probVaccWorks, excessCFR, col = factor(cfr), group=cfr))  + 
        geom_hline(yintercept=0, linetype = 2, col = 'black') + guides(col = guide_legend(rev=T, title='CFR') ) +
            ylab('excess risk of death in control arm') + xlab('probability vaccine works')  + 
                ylim(-.01, .03) + scale_x_continuous(breaks=seq(0,1,b=.2))
print(p)
graphics.off()

ps <- 1/200
tmp <- eqd[probSAE==ps & infRisk==.05]
p <- ggplot(tmp) + facet_grid(~efficacy) + scale_color_manual(values=cols) + # palette='RdBu') + 
    geom_line(aes(probVaccWorks, excessCFR, col = factor(cfr), group=cfr))  + 
        geom_hline(yintercept=0, linetype = 2, col = 'black') + guides(col = guide_legend(rev=T, title='CFR')) +
            ylab('excess risk of death in control arm') + xlab('probability vaccine works')  + 
                ggtitle(paste0('pSAE = ',ps)) + scale_y_continuous(labels=percent)#, limits=c(-.01,.06))
print(p)
