if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls4', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('nid', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/EbolaVaccSim/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')
rm(list=ls(all=T)); gc()
require(optiRum); require(RColorBrewer); require(boot); require(data.table); require(ggplot2); require(grid); require(reshape2); require(parallel)
## Multiple individual priors
setwd('Figures')
newdr <- 'Fig1Uncertainty'
if(!dir.exists(newdr)) dir.create(newdr)
setwd(newdr)

shadeplot <- function(x,y, ...) {
polygon(c(x,rev(x)), c(y,rep(0,length(y))), ...)
}


effPrior <- function(efficacyMean=.75, probEff = .5, xlim = c(.05,20), 
                     xtcks = c(.05,.1,.25,.5,1,2,4,10,20), ylim = c(0,1.2),
                     scale1 = .6, scale2 = .2,  xlab = '', par = list(mar=rep(0,4)),
                     showAxes=T, showav=T, main = '', saveNm=NA, w=1, h=.7)
{
    x <- seq(log(xlim[1]), log(xlim[2]), l = 200)
    y <- probEff*dcauchy(x, location =log(1-efficacyMean), scale=scale1) + (1-probEff)*dcauchy(x, location=log(1), scale=scale2)
    xats <- log(xtcks)
    if(!is.na(saveNm)) {
        pdf(paste0(saveNm,'.pdf'), w=h, h=h)
        par(par)
    }
    plot(x,y, type='n', xlim = range(x), axes=F, ylab='', xlab =xlab, ylim = ylim, main =main)
    shadeplot(x,y, col = 'black')
    margMean <- probEff*log(1-efficacyMean) + (1-probEff)*log(1)
    if(showav) segments(margMean,0,margMean,ylim[2]/3, lwd = 3, col = 'red')
    if(showAxes) {
        axis(1, at = xats, labels=xtcks)
        axis(1, at = xats[xtcks<=1], labels=paste0((1-xtcks[xtcks<=1])*100,'%'), line = 2, lty = 0)
    }
    if(!is.na(saveNm)) dev.off() 
}



saePrior <- function(saeMean=-4, scale = 1.3,
                     xtcks = -6:-2, ylim = c(0, 1.2), calcYmax = F,
                     scale1 = .6, xlab = '', browse = F,
                     showAxes=T, main = '', saveNm=NA, w=1, h=.7)
{
    if(browse) browser()
    xlim <- range(xtcks)
    x <- seq(min(xtcks), max(xtcks), l = 200)
    y <- dcauchy(x, location =saeMean, scale=scale1)
    #xats <- paste0('10^', xtcks)
    xats=parse(text=paste("10", "^", xtcks))
    if(!is.na(saveNm)) {
        pdf(paste0(saveNm,'.pdf'), w=h, h=h)
     #   par(mar=rep(0,4))
    }
    if(calcYmax) ylim[2] <- max(y)
    plot(x,y, type='n', xlim = range(x), axes=F, ylab='', xlab =xlab, ylim = ylim, main =main)
    shadeplot(x,y, col = 'black')
    if(showAxes) {
        axis(1, at = xtcks, labels=xats)
    }
    if(!is.na(saveNm)) dev.off() 
}

pdf('sae.pdf', w = 3, h = 2)
par(mar =c(4,.5,1,1), 'ps'=16)
saePrior(saeMean=-4.5, scale1=.6, xlab = 'probability of SAE', calcYmax=T)#, saveNm='sae')
dev.off()
     
pdf(file='3 panel.pdf', w = 5.5, h = 2)
par(bty ='n', mar = c(7,.5,2,0), mfrow = c(1,3))
effPrior(probEff=0, scale2 = 1, main = 'agnosticism') ## pure agnosticism on efficacy
effPrior(probEff=.7,main= 'optimistic conflict') ## strong conflict on efficacy
effPrior(probEff=.3, scale1 = .8, main='pessimistic conflict') ## weak conflict on efficacy
dev.off()

effPrior(probEff=.2, efficacyMean = .8, scale1= .3, scale2 = .8, showav = F, showAxes = F, saveNm='1')
effPrior(probEff=.7, scale2 = .3,  showav = F, showAxes = F, saveNm='2')
effPrior(probEff=.3, scale1 = .8, scale2 = .3, showav = F, showAxes = F, saveNm='3')
effPrior(probEff=.1, efficacyMean = .7, scale1= .4, scale2 = .8, showav = F, showAxes = F, saveNm='4') 
effPrior(probEff=.6, efficacyMean = .8, scale1= .2, scale2 = .8, showav = F, showAxes = F, saveNm='5') 
effPrior(probEff=.2, efficacyMean = .8, scale1= .4, scale2 = .9, showav = F, showAxes = F, saveNm='6')

pdf(file='3 panel.pdf', w = 5.5, h = 2)
par(bty ='n', mar = c(7,.5,2,0), mfrow = c(1,3))
effPrior(probEff=0, scale2 = 1, main = 'agnosticism') ## pure agnosticism on efficacy
effPrior(probEff=.2, efficacyMean = .8, scale1= .4, scale2 = .9, showav = F, showAxes = T)
dev.off()


pdf(file='expertComm.pdf', w = 3, h = 2)
par(bty ='n', mar = c(5,0,0,0))
effPrior(probEff=.2, efficacyMean = .8, scale1= .4, scale2 = .9, showav = F, showAxes = T)
graphics.off()
## SAEs
  
