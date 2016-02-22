set.seed(1)
rfs <- rlnorm(100)

simTrial <- function(s=1, het = F) {
    dt <- data.table(expand.grid(id=1:100, week = 1:50))
    ## dt <- merge(dt, data.table(id=1:100, haz = rlnorm(100)), by = 'id')
    if(het) indivRR <- rfs else indivRR <- 1
    indivTab <- data.table(id=1:100, haz = .03, indivRR= indivRR, iHaz=.03*indivRR)
    dt <- merge(dt, indivTab, by = 'id')
    dt$s <- s ## simulation #
    dt$infected <- 0L
    naive <- T ## has not had end criteria satisfied
    for(ii in 1:50) {
        infected <- dt[infected>0, unique(id)] ## whose already been infected
        uninfected <- (1:100)[!1:100 %in% infected] ## who hasn't 
        dt[id %in% uninfected & week==ii, infected:=rbinom(.N, 1, 1-exp(-iHaz))] ## run infection for those who haven't
        if(naive & dt[infected>0, length(unique(id))] > 30) { ## once hitting 30 infections reduce hazard by 90% (only 1st time, i.e. naive==T)
            dt[week > ii, iHaz:= iHaz * 0.1] ## reduce individual hazards by 90% 
            naive <- FALSE
        }
    }
    dt[infected>0, length(unique(id))/100]
    return(dt)
}

simNtrials <- function(nn, het) {
    blah <- list()
    length(blah) <- nn
    for(ss in 1:nn) {
        blah[[ss]] <- simTrial(ss, het = het)
    }
    blah <- rbindlist(blah)
    setkey(blah, s, id)
    return(blah)
}

## Pure simulation based approach
collectRisks <- function(blah) {
    infs <- blah[,list(inf = as.numeric(sum(infected)>0)),list(id,s)]
    simRisk <- infs[,list('risk'=mean(inf)),id]
    return(simRisk)
}

collectHazards <- function(blah) {
    ## **DOES NOT WORK** right-censoring: for weeks after someone has been infected, set wt to 0
    ##    blah$wt <- 1
    ##    blah[,wt:=1-as.numeric(1:50>which(infected==1)),list(id,s)]


    blah[,hazC:=iHaz*wt] ## create censored hazards based on censorship
    ## Option 1: average hazards, then calculate risk
    simHazH <- blah[,list('haz'=sum(hazC)),list(id,s)] ## sum censored hazards for each individual-simulation
    simHazH <- simHazH[, list(haz=mean(haz)), id] ## take mean of summed censored hazards for individuals across simulations
    simHazH[,riskH_1:=1-exp(-haz)]                 ## convert those mean summed censored hazards to risks
    ## Option 2: calculate risk, then average risks
    simHazR <- blah[,list('risk'=1-exp(-sum(hazC))),list(id,s)] ## sum censored hazards for each individual-simulation & convert to risk
    simHazR <- simHazR[, list(riskH_2=mean(risk)), id] ## take mean of summed censored hazards for individuals across simulations
    simHaz <- merge(simHazH, simHazR, by = 'id')
    return(simHaz)
}

## homogenous version
blah <- simNtrials(1000, F)
cR <- collectRisks(blah)
cR[,mean(risk)] ## mean risk averaged across homogenous individuals
cH <- collectHazards(blah)
comp <- merge(cR, cH, by = 'id')
comp

## heterogeneous version
blahHet <- simNtrials(1000, T)
cR <- collectRisks(blahHet)
## cR[,mean(risk)] ## mean risk averaged across homogenous individuals (not relevant for heterogeneous sim)
cH <- collectHazards(blahHet)
compHet <- merge(cR, cH, by = 'id')
compHet

pdf('Figures/testJD.pdf')
with(compHet, plot(risk, riskH_1))
rp <- function(x) 1-exp(-x)
yx <- function(x) x
curve(rp, add=T)
curve(yx, add=T)
graphics.off()
