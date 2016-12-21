
percent <- function(x) paste0(formatC(x*100), '%')

gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}
group.colors <- c(RCT = "#333BFF", FRCT = "#CC6600", SWT ="#9633FF")
group.colors[c(1,3,2)] <- gg_color_hue(3)
group.colors['SWCT'] <- 'orange'
if(class(pf)!='function') {
pf$trial <- factor(pf$trial, levels=levels(pf$trial)[c(2,1,3)])
pf[, biasNAR:=biasNAR/vaccEff]
levels(pf$order)[1] <- 'random'
levels(pf$order)[2] <- 'risk-prioritized'
}

####################################################################################################
## ggplot theme for MS
thax <- element_text(colour = 'black', size = 8)
thsb <- theme(axis.text.x = thax, axis.text.y = thax, plot.title = element_text(vjust=1),
              axis.title.y = element_text(vjust = 1), axis.title.x = element_text(vjust = -.5),
              axis.line = element_line(), axis.ticks = element_line(color='black'),
              panel.spacing = unit(1, "lines"), legend.key.height=unit(1.3,"line")
              , strip.background = element_rect(fill = NA)
              ,legend.position = 'right'
              #, axis.line = element_blank()
              ,panel.grid.major = element_blank()
              , panel.grid.minor = element_blank()
              ,panel.border = element_blank()
              ,panel.background = element_blank()
            , legend.background =  element_blank()
            , legend.key =  element_blank()
            , legend.key.width=unit(2.5,"line")
            #,legend.justification=c(1,1), legend.position=c(1,1)
            # ,legend.position='top'
              )
theme_set(theme_grey(base_size = 12))

eb <- theme(
    ## axis.line=element_blank(),axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),legend.position="none",
    panel.background=element_blank(),
    panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),plot.background=element_blank())


makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
