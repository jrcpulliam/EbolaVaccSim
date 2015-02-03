## Check end trial functions work. Complicated options and ordering of who gets vaccinated in what order.

## ## none
## p1 <- simTrial(makeParms('CRCT',small=F, ord='none'))
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## t1 <- seqStop(s1, verbose=0)
## t1 <- endCRCT(t1)
## print(t1$endTrialDay)
## with(t1, ## vaccDay should increment as time-updated clusHaz (order given in V5)
##      print(EVpopH[cluster %in% unVaccClusters &idByClus==1 & day == endTrialDay,
##                   list(cluster, vaccDay, day, clusHaz, order(rev(order(clusHaz))))])
##      )

## ## BL
## p1 <- simTrial(makeParms('CRCT',small=F, ord='BL'))
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## t1 <- seqStop(s1, verbose=0)
## t1 <- endCRCT(t1)
## print(t1$endTrialDay)
## with(t1, ## vaccDay should increment as time-updated clusHaz (order given in V5)
##      print(EVpopH[cluster %in% unVaccClusters &idByClus==1 & day == endTrialDay,
##                   list(cluster, vaccDay, day, clusHaz, order(rev(order(clusHaz))))])
##      )

## ## TU
## p1 <- simTrial(makeParms('CRCT',small=F, ord='TU'))
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## t1 <- seqStop(s1, verbose=0)
## t1 <- endCRCT(t1)
## print(t1$endTrialDay)
## with(t1, ## vaccDay should increment as time-updated clusHaz (order given in V5)
##      for(ii in 1:length(unVaccClusters)) 
##      print(EVpopH[cluster %in% unVaccClusters &idByClus==1 &day ==vaccDaysLeft[ii], 
##                   list(cluster, vaccDay, day, clusHaz, order(rev(order(clusHaz))))])
##      )


## ## none, option 1 
## p1 <- simTrial(makeParms('RCT',small=F, ord='none', RCTendOption = 1))
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## t1 <- seqStop(s1, verbose=0)
## t1 <- endRCT(t1)
## with(t1, 
##      print(EVpopH[idByClus %in% c(1,clusSize/2+1) & day == endTrialDay,
##                   list(cluster, vaccDay, endTrialDay = day, clusHaz, ord = order(rev(order(clusHaz))))])
##      )

## ## none, option 2
## p1 <- simTrial(makeParms('RCT',small=F, ord='none', RCTendOption = 2))
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## t1 <- seqStop(s1, verbose=0)
## t1 <- endRCT(t1)
## with(t1, 
##      print(EVpopH[idByClus %in% c(1,clusSize/2+1) & day == endTrialDay,
##                   list(cluster, vaccDay, endTrialDay = day, clusHaz, ord = order(rev(order(clusHaz))))])
##      )

## ## BL, option 1
## p1 <- simTrial(makeParms('RCT',small=F, ord='BL', RCTendOption = 1))
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## t1 <- seqStop(s1, verbose=0)
## t1 <- endRCT(t1)
## with(t1, EVpopH[idByClus == (clusSize/2+1) & day == endTrialDay, ## vacc group w/in each cluster
##           list(cluster, vaccDay, endTrialDay = day, clusHaz, ord = order(rev(order(clusHaz))))])
## with(t1, EVpopH[idByClus %in% 1 & day == endTrialDay,  ## control group w/in each cluster
##           list(cluster, vaccDay, endTrialDay = day, clusHaz, ord = order(rev(order(clusHaz))))])

## ## BL, option 2
## p1 <- simTrial(makeParms('RCT',small=F, ord='BL', RCTendOption = 2))
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## t1 <- seqStop(s1, verbose=0)
## t1 <- endRCT(t1)
## with(t1, EVpopH[idByClus == (clusSize/2+1) & day == endTrialDay, ## vacc group w/in each cluster
##           list(cluster, vaccDay, endTrialDay = day, clusHaz, ord = order(rev(order(clusHaz))))])
## with(t1, EVpopH[idByClus %in% 1 & day == endTrialDay,  ## control group w/in each cluster
##           list(cluster, vaccDay, endTrialDay = day, clusHaz, ord = order(rev(order(clusHaz))))])


## ## TU, option 1
## p1 <- simTrial(makeParms('RCT',small=F, ord='TU', RCTendOption = 1))
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## t1 <- seqStop(s1, verbose=0)
## t1 <- endRCT(t1)
## with(t1, 
##      print(EVpopH[idByClus %in% (clusSize/2+1) & day == endTrialDay, ## vacc group w/in each cluster
##             list(cluster, vaccDay, endTrialDay = day, clusHaz, ord = order(rev(order(clusHaz))))])
##      )
## t1$endTrialDay
## with(t1, 
##      for(ii in 1:length(notYetVaccClusters)) 
##      print(EVpopH[idByClus %in% 1 & day == endTrialDay + (ii-1)*delayUnit,  ## control group w/in each cluster
##             list(cluster, vaccDay, day, clusHaz, ord = order(rev(order(clusHaz))))])
##      )

## ## TU, option 2
## p1 <- simTrial(makeParms('RCT',small=F, ord='TU', RCTendOption = 2))
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## t1 <- seqStop(s1, verbose=0)
## t1 <- endRCT(t1)

## with(t1, 
##      #for(ii in 1:length(notYetVaccClusters)) 
##      print(EVpopH[idByClus %in% (clusSize/2+1) & day == endTrialDay,# + (ii-1)*delayUnit, ## vacc group w/in each cluster
##             list(cluster, vaccDay, endTrialDay = day, clusHaz, ord = order(rev(order(clusHaz))))])
##      )

## t1$endTrialDay
## with(t1, 
##      for(ii in 1:length(notYetVaccClusters)) 
##      print(EVpopH[idByClus %in% 1 & day == endTrialDay + (ii-1)*delayUnit,  ## control group w/in each cluster
##             list(cluster, vaccDay, day, clusHaz, ord = order(rev(order(clusHaz))))])
##      )

####################################################################################################
## Look at trials without any logistical delays
## p1 <- simTrial(makeParms('RCT',small=F, ord='TU', RCTendOption = 2, delayUnit = 0))
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## t1 <- seqStop(s1, verbose=0)
## t1 <- endT(t1)

## with(t1, 
##      print(EVpopH[idByClus %in% c(1,clusSize/2+1) & day == endTrialDay,
##                   list(cluster, vaccDay, endTrialDay = day, clusHaz, ord = order(rev(order(clusHaz))))])
##      )

## p1 <- simTrial(makeParms('CRCT',small=F, ord='TU', RCTendOption = 2, delayUnit = 0))
## s1 <- makeSurvDat(p1)
## s1 <- activeFXN(s1)
## t1 <- seqStop(s1, verbose=0)
## t1 <- endT(t1)

## with(t1, ## vaccDay should increment as time-updated clusHaz (order given in V5)
##      print(EVpopH[idByClus==1 & day == endTrialDay,
##                   list(cluster, vaccDay, day, clusHaz, order(rev(order(clusHaz))))])
##      )
