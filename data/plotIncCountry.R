if(grepl('Stevens-MBP', Sys.info()['nodename'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('stevebellan', Sys.info()['login'])) setwd('~/Documents/R Repos/EbolaVaccSim/')
if(grepl('ls', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/VaccEbola/')
if(grepl('wrang', Sys.info()['nodename'])) setwd('/home/02413/sbellan/work/EbolaVaccSim/')

require(gdata); require(data.table); library(RColorBrewer); library(mgcv);  require(RColorBrewer)
ff <- list.files('data/Data',full.names = T)

for(ii in 1:3) assign(paste0('dat',ii), as.data.table(read.csv(ff[ii])))

dat <- rbind(dat1,dat2,dat3)

summary(dat)

pdat <- dat[Ebola.measure=='Number of cases' & Ebola.data.source=='Patient database']


pdat[, dateWk:=sub(' \\(.+)', '', sub('.+to ', '', Epi.week))]
pdat[, date:=as.Date(dateWk, format='%d %B %Y')]
pdat <- pdat[order(Location, date)]

ebinc <- function(upto='2015-11-01', bg='black', fg='white',ps=27, detail=T, cex = 1) {
    upto <- as.Date(upto)
    ## plot all 3 countries subnational
    countries <- pdat[,unique(Country)]
    par(mfrow=c(3,1), mar = c(6,5,1,0), oma = c(3.5,10,0,0), ps = ps, bg = bg, col.axis=fg, col.sub = fg, col.lab = fg, col.main = fg)
    for(ci in 1:3) { 
        cc <- countries[ci]
        tempmax <- pdat[Location!='' & Case.definition=='Confirmed'  & !is.na(Numeric) & Country==cc, max(Numeric)]
        pdat[, plot(date,Numeric, , type = 'n', bty = 'n', axes=F, ylab='', xlab='', las = 1, ylim = c(0,tempmax))]
        ## wktcks <- seq.Date(min(pdat[,date],na.rm=T),max(pdat[,date],na.rm=T), by='week')
        wktcks <- seq.Date(as.Date('2013-12-01'), as.Date('2015-11-01'), by='month')
        wklabs <- wktcks
        wklabs[wklabs>upto] <- NA
        ## wklabs[1:length(wklabs)%%4!=0] <- NA
        axis.Date(1, at=wktcks, las = 2, lab = '', col=fg)
            axis(2, pretty(c(0,tempmax)), col=fg, las = 1)
            axis.Date(1, at=wktcks, las = 2, lab = format(wklabs, format='%b'), cex=2, col=fg)
        country.col <- data.table(loc= pdat[Location!='' & Country==cc, unique(Location)])
        ## country.col[, col:=rainbow(nrow(country.col))]
        country.col[, col:=brewer.pal(12, 'Set3')[1:nrow(country.col)%%9 + 1]]
        pdat[date <= upto & Location!='' & Case.definition=='Confirmed'  & !is.na(Numeric) & Country==cc, ## & Location=='' 
             lines(date, Numeric, col = country.col[loc==Location,col], lwd = 2), Location]
        abline(h=0, col = 'dark gray', lty = 2)
        ##    legend('topright', leg = country.col$loc, col = country.col$col, lwd = 1, bty = 'n', ncol=2)
        pdat[!is.na(Numeric) & Country==cc& Location!='', list(date, Location, Numeric,Ebola.measure), Location]
        mtext(cc, side = 3, -6, adj=.1, col=fg, cex = cex)
        if(detail) {
            if(ci==3) {
                axis.Date(1, c(as.Date('2014-01-01'),as.Date('2014-12-31')), line=5)
                axis.Date(1, c(as.Date('2015-01-01'),as.Date('2014-12-31')), line=5)
            }
            ## Liberia
            if(cc=='Liberia' & upto>=as.Date('2015-02-02')) {
                rect(as.Date('2015-02-02'), 0, min(upto, as.Date('2015-09-02')), 60, col=gray(.9,.5))
                text(as.Date('2015-02-02'), 75, "NIH vaccine trial", col=fg, pos = 4)
            }
            if(cc=='Liberia' & upto>=as.Date('2015-09-02')) {
       #         arrows(as.Date('2015-09-02'), 35, as.Date('2015-09-02'), 5, len = .05, col='yellow', lwd = 3)
       #         text(as.Date('2015-09-02'), 35, "halted", col='yellow', pos = 4)
            }
            ## SL
            if(cc=='Sierra Leone' & upto>=as.Date('2015-04-15')) {
                rect(as.Date('2015-04-15'), 0, min(upto, as.Date('2015-08-31')), 60, col=gray(.9,.5))
                text(as.Date('2015-04-15'), 75, "CDC vaccine trial", col=fg, pos = 4)
            }
            if(cc=='Sierra Leone' & upto>=as.Date('2015-08-31')) {
          #      arrows(as.Date('2015-08-31'), 35, as.Date('2015-08-31'), 5, len = .05, col='yellow', lwd = 3)
           #     text(as.Date('2015-08-31'), 35, "halted", col='yellow', pos = 4)
            }
            ## Guinea
            if(cc=='Guinea' & upto>=as.Date('2015-03-25')) {
                rect(as.Date('2015-03-25'), 0, min(upto, as.Date('2015-07-31')), 25, col=gray(.9,.5))
                text(as.Date('2015-03-25'), 45, "WHO vaccine trial", col=fg, pos = 4)
            }
            if(cc=='Guinea' & upto>=as.Date('2015-07-31')) {
                arrows(as.Date('2015-07-31'), 25, as.Date('2015-07-31'), 5, len = .05, col='yellow', lwd = 3)
                text(as.Date('2015-07-31'), 25, "found\nefficacious", col='yellow', pos = 4)
            }
            ## vacc plans
            par(xpd=NA)
            if(ci==2  & upto>=as.Date('2014-09-30')) text(as.Date('2014-09-30'), 180, 'vaccine\ndiscussions\nbegin', col='yellow')
            par(xpd=F)
        }
    }
    mtext(expression(frac(cases,week)), side=2, line=0.3, outer=T, col=fg, las = 1, cex = 1)
}


## ebinc(wktcks[1], ps = 16)

wktcks <- as.character(seq.Date(min(pdat[,date],na.rm=T),max(pdat[,date],na.rm=T), by='week'))
## wktcks <- wktcks[80:85]

## Movie
require(animation);
resScl <- 1.5
nm <- paste0('ebolaIncnodetail.mov')
if(file.exists(nm)) file.remove(nm)
saveVideo({
    ani.options(interval = 0.15, nmax = 300, ani.dev='png', ani.type='png')
    for(ww in wktcks) ebinc(ww, ps=30, detail=F )
}, video.name = nm, other.opts = "-b 3000k -pix_fmt yuv420p", ani.width = 800*resScl, ani.height = 600*resScl)

## resScl <- 1
## png('EbolaIncXDistrictWA Tall.png', w=800*resScl, h=1000*resScl)
 ebinc(ps = 30, detail=F, cex = 1)
## dev.off()

evddat <- pdat[!is.na(Numeric) & Case.definition=='Confirmed' & Location!='',.(date, Numeric, Country,Location)]
setnames(evddat, 'Numeric','cases')

SLdistr <- evddat[Country=='Sierra Leone']
SLcountry <- SLdistr[,.(cases = sum(cases)), date]
Gdistr <- evddat[Country=='Guinea']
Gcountry <- Gdistr[,.(cases=sum(cases)), date]
Ldistr <- evddat[Country=='Liberia']
Lcountry <- Ldistr[,.(cases=sum(cases)), date]

nms <- c('evddat','SLdstr','SLcountry','Gdistr','Gcountry','Ldistr','Lcountry')
for(nn in nms) assign(nn, data.frame(get(nms)))

save(evddat, SLdistr, SLcountry, Gdistr, Gcountry, Ldistr, Lcountry, guecDat, file='data/WAevddat.Rdata')

pdf('data/All EVD.pdf', w=10, h = 8)
par(ps=16, bty='n')
cols <- c('red','purple','orange')
evddat[, plot(date,cases, type = 'n', ylab = 'weekly EVD incidence by district', xlab='', xaxt='n')]
wktcks <- seq.Date(as.Date('2013-12-01'), as.Date('2015-11-01'), by='month')
axis.Date(1, at=wktcks, las = 2, lab = format(wktcks, format='%b-%y'), cex=2)
evddat[,lines(date,cases, col = cols[Country], lwd=2), Location]
legend('topright', leg = evddat[,levels(Country)], col=cols, pch = 16)
graphics.off()


for(cc in levels(evddat$Country)) {
    pdf(paste0('data/', cc, ' EVD.pdf'), w=10, h = 8)
    par(ps=16, bty='n')
    temp <- evddat[Country==cc]
    temp[,Location:= levels(Location)[as.numeric(Location)]]
    temp[,Location:=factor(Location)]
    nLoc <- temp[,length(unique(Location))]
    temp[, colx:= rainbow(nLoc)[as.numeric(Location)]]
    evddat[, plot(date,cases, type = 'n', ylab = 'weekly EVD incidence by district', xlab='', xaxt='n')]
    wktcks <- seq.Date(as.Date('2013-12-01'), as.Date('2015-11-01'), by='month')
    axis.Date(1, at=wktcks, las = 2, lab = format(wktcks, format='%b-%y'), cex=2)
    for(ll in levels(temp[,Location])) temp[Location==ll, lines(date,cases, col = colx, lwd=2)]
    temp[!duplicated(Location), legend('topright', leg = Location, col=colx, pch = 16, ncol = 2)]
    graphics.off()
    pdf(paste0('data/', cc, ' EVD district per page.pdf'), w=10, h = 8)
    for(ll in levels(temp[,Location])) {
    evddat[, plot(date,cases, type = 'n', ylab = 'weekly EVD incidence by district', xlab='', xaxt='n', main = ll)]
    axis.Date(1, at=wktcks, las = 2, lab = format(wktcks, format='%b-%y'), cex=2)
    temp[Location==ll, lines(date,cases, lwd=2)]
    }
    graphics.off()
}
