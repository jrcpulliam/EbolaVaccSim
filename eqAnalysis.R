if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('nid', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/EbolaVaccSim/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')

####################################################################################################
## Analytical equipoise functions
require(data.table); require(ggplot2); require(scales); require(gsDesign); require(pwr); require(RColorBrewer); require(grid)
source('ggplotTheme.R')
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

####################################################################################################
## excess risk, size, infection risk, power
dat <- CJ(infRisk = seq(.0001, .2, l = 100), size = seq(10, 10^4, by = 10),
          probSAE = 10^-4, efficacy=.8, cfr = .7, probVaccWorks = .5)

dat <- CJ(loginfRisk = seq(log(.0001), log(.2), l = 100), logsize = seq(log(10), log(10^4), l = 100),
          probSAE = 10^-4, efficacy=.8, cfr = .7, probVaccWorks = .5)
dat[,size:=exp(logsize)]
dat[,infRisk:=exp(loginfRisk)]

dat[,power:=power.prop.test(size, infRisk, infRisk*(1-efficacy))$power]
dat[, vaccCFR := probSAE + (1-probSAE) * (probVaccWorks * (infRisk * (1-efficacy) * cfr) + (1-probVaccWorks) * infRisk * cfr)]
dat[, contCFR := infRisk * cfr]
dat[, excessCFR := contCFR - vaccCFR]

v <- ggplot(dat, aes(x=infRisk, y=size, z=power))
v + geom_raster(aes(fill = excessCFR)) + thsb + 
  geom_contour(aes(z= power,colour = ..level..), size = 2) +
      scale_y_continuous(breaks=c(10,30,100,300,1000,3000,10000), trans="log10") +
          scale_x_continuous(breaks=c(.001, .003, .01, .03, .1, .3), trans="log10") +
      scale_fill_continuous(low='blue', high="red") +
      scale_color_gradient(low='orange', high="white")          +
          xlab('baseline infection risk') + ylab('trial power')
ggsave(file.path('Figures','size risk power.pdf'), w = 6, h = 4)

scl <- 10
dat <- expand.grid(x = scl * seq(0, 1, by = 0.01), 
                   y = scl * seq(0, 1, by = 0.01))
dat$z <- ((scl - dat$x) * dat$y) / ((scl - dat$x) * dat$y + 1)

# create the plot, the geom_contour may not be needed, but I find it helpful
ggplot(dat) + 
aes(x = infRisk, y = size, z = power, fill = power) + 
geom_tile() + 
geom_contour(color = "white", alpha = 0.5) 
