## census data 2004
# district N
# n.b. sherbro urban merged into Bonthe
load(file='data/cleanSLData.Rdata')
pop.dt <- data.table(reg = factor(levels(sl$reg)), N=c(
 463668,
 408390,
 139687,
 358190,
 270462,
 497948,
 265758,
 335401,
 260910,
 453746,
 228392,
 347197,
 174249,
 772873
), key="reg")

pop.dt[,proportion := N/sum(N)]
setkey(sl,"Date","reg")
pop.dt[sl[,list(total_inc = sum(cases)), by="reg"], cumulative_inc := total_inc ]
pop.dt[, inc_per_capita := cumulative_inc / N ]

ggplot(pop.dt) + theme_bw() + aes(x=N, y=decays, color=reg, group=1) + geom_point() + stat_smooth(geom = "line", method = "glm")

## Fit to current incidence trends
source('ExpFit.R');
fits <- NULL
for(rr in levels(sl$reg)) fits[[rr]] <- doProj(sl[reg==rr], ll='exp_nbinom_ll')
fit.dt <- data.table(reg=factor(names(fits)), decay_rates = sapply(fits,function(fit) fit$fit$par['decay_rate']))

pop.dt[fit.dt, decays := decay_rates]

samp.size <- 10000

test <- sample(pop.dt$decays, replace = T, size = samp.size)
test.prop <- sample(pop.dt$decays, replace = T, size = samp.size, prob = pop.dt$proportion)

sd(test)
sd(test.prop)

test.dt <- data.table(type=factor(rep(c("unif","prop"), each=samp.size)), samp=c(test, test.prop) )

ggplot(test.dt) + aes(group=type, fill=type, x=samp) + geom_bar(position = "dodge")
