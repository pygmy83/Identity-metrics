# FigMutiMetrics
# library(ggplot2), library(gridExtra) 

# https://stats.stackexchange.com/questions/177129/ggplot-and-loops
# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html


oneplotwidth <- 1.4
oneplotheight <- 1.4

varstoplotX <- c(1:5)
#varstoplotXnames <- names(dataz)[1:5]
varstoplotXnames <- c('Individuality', 'Variables', 'Covariance', 'Individuals', 'Observations')
varstoplotY <- c(7:10)
varstoplotYnames <- names(dataz)[7:10]
xypairs <- data.frame(param=rep(varstoplotX, length(varstoplotY)),
                      metric=rep(varstoplotY, each=length(varstoplotX)))
xypairs.names <- data.frame(param=rep(varstoplotXnames, length(varstoplotY)),
                            metric=rep(varstoplotYnames, each=length(varstoplotX)))
plots <- list()
results <- matrix(ncol=8, nrow=nrow(xypairs.names))

for (i in 1:nrow(xypairs)) {
  param <- as.numeric(xypairs[i,1])
  metric <- as.numeric(xypairs[i,2])
  dataz.temp <- data.frame(x=as.factor(dataz[,param]), y=dataz[,metric])
  
  # get numerical results of comaprisons
  result.aov <- aov(y ~ as.factor(x), data=dataz.temp)
  result.Tukey <- data.frame(TukeyHSD(result.aov)$'as.factor(x)')
  summary(result.aov)
  results[i,1] <- paste(xypairs.names[i,2], "X", xypairs.names[i,1])
  results[i,2] <- summary(result.aov)[[1]][1,1]
  results[i,3] <- summary(result.aov)[[1]][2,1]
  results[i,4] <- round(summary(result.aov)[[1]][1,4], 2)
  results[i,5] <- round(summary(result.aov)[[1]][1,5], 2)
  
  wanted.n <- nlevels(dataz.temp$x)-1
  wanted <- rep(NA, wanted.n)
  wanted[1] <- 1
  for (j in 2:wanted.n){
    wanted[j] <- wanted[j-1]+wanted.n
    wanted.n <- wanted.n-1
  }
    
  comparison <- row.names(result.Tukey)[wanted]
  pvals <- round(result.Tukey$p.adj[wanted],3)
  results[i,6] <- ifelse(length(comparison[pvals>0.05])==0, "all sig", paste(" ", paste(comparison[pvals<0.05],collapse="  ")))
  results[i,7] <- ifelse(length(tail(comparison[pvals<0.05],1))==0, 'no sig', paste(" ", comparison[pvals<0.05][1]))
  results[i,8] <- ifelse(length(tail(comparison[pvals<0.05],1))==0, 'no sig', paste(" ", tail(comparison[pvals<0.05],1)))
  
  dataz.temp <- data.frame(x=dataz[,param], y=dataz[,metric])
  dataz.sum <- summarySE(dataz.temp, measurevar='y', groupvars='x')
  #xvalues <- theme(axis.text.x = element_text(color="black"))
  # xvalues <- theme(axis.text.x = element_blank())
  
  # whether to plot x lab
  if (i %in% c(16, 17, 18, 19,20)) {
    #xtitle <- theme(axis.title.x = element_text(color="black"))
    xlabname <- as.character(xypairs.names[i,1])
  } else {
    xlabname <- as.character(xypairs.names[i,1]) #''
  }
  
  # whether to plot y lab
  if (i %in% c(1, 6, 11, 16)) {
    #xtitle <- theme(axis.title.x = element_text(color="black"))
    ylabname <- as.character(xypairs.names[i,2])
  } else {
    ylabname <- as.character(xypairs.names[i,2]) #''
  }
  
  
  plots[[i]] <- 
    ggplot(dataz.sum, aes(x=x, y=y)) + 
    geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0) +
    geom_line() +
    geom_point(size=1, shape=21, fill="black") +
    theme_bw() + 
    xvalues + 
    theme(text = element_text(size=8, color='black'), 
          #axis.title.x = element_blank(), 
          #axis.title.y = element_blank(),
          #axis.text.x = element_text(color="black"), #element_blank(), # 
          #axis.text.y = element_blank(), #element_text(color="black"), #  
          #axis.ticks = element_line(size = 2),
          axis.ticks.length=unit(.10, "cm"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()
    ) +  
    xlab(xlabname) + # xlab(names(dataz)[param]) +
    ylab(ylabname) + # ylab(names(dataz)[metric]) +
    ggtitle("")  
  #expand_limits(y=c(dataz.sum$min, dataz.sum$max))
  #scale_y_continuous(breaks=0:4*0.25)
}

results <- data.frame(results)
names(results) <- c('pair', 'df', 'dfres', 'F', 'p', 'sigTukeys', 'firstsigTukey', 'lastsigTukey')
write.table(results, 'res-IterateDatasetsMulti.csv', sep=';',row.names=F)

#margin = theme(plot.margin = unit(c(2,2,2,2), "cm"))
#grid.arrange(grobs = lapply(pl, "+", margin))

win.metafile("mygraph.wmf", width=length(varstoplotX)*oneplotwidth, height=length(varstoplotY)*oneplotwidth)
grid.arrange(grobs=plots, ncol=length(varstoplotX))
dev.off()



### ANOTHER VARIANT WITH GPLOTS:

# plot figures for Fig XX 'Multivariate identity metrics'
varstoplotX <- c(1:5)
varstoplotXnames <- names(dataz)[1:5]
varstoplotY <- c(7:10)
varstoplotYnames <- names(dataz)[7:10]

pdf("mygraph.pdf", width=7, height=6)
par(mfrow=c(4,5), mar=c(0,0,0,0), oma=c(3,3,0.1,0.1))
for (i in varstoplotY) {
  for (j in varstoplotX) {
    y <- dataz[,i]
    x <- dataz[,j]
    plotmeans(y~x, bars=T,n.label=F,barcol='black',ylim=c(min(y), max(y)),
              yaxt ='n', xaxt='n',
              cex=1, cex.lab=0.8, cex.axis=0.8)
    #xlab=names(dataz)[j],
    #ylab=names(dataz)[i], 
    #    summary(aov(y~as.factor(x), data=dataz))
    #    TukeyHSD(aov(y~as.factor(x), data=dataz))
  }
}
dev.off()

