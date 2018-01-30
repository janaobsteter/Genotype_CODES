#selection eff po replikah grafi



#to je max-min po ponovitvah, zgoraj je povprečje
maxmin <- data.frame(scenario=NA, rep = NA, minGenicSD=NA, maxGenicSD=NA, minTGV=NA, maxTGV=NA)
row = 1
for (rep in unique(TGVsAll$Rep)) {
  for (scenario in unique(TGVsAll$scenario)) {
    maxmin[row,] <- c(scenario, rep, min(TGVsAll$SDGenicSt[TGVsAll$scenario==scenario & TGVsAll$Rep == rep]), max(TGVsAll$SDGenicSt[TGVsAll$scenario==scenario & TGVsAll$Rep == rep]), min(TGVsAll$zMeanGenic[TGVsAll$scenario==scenario & TGVsAll$Rep == rep]), max(TGVsAll$zMeanGenic[TGVsAll$scenario==scenario & TGVsAll$Rep == rep]))
    row <- row +1
  }
}


#genska
maxmin$minGenicSD <- as.numeric(maxmin$minGenicSD)
maxmin$maxGenicSD <- as.numeric(maxmin$maxGenicSD)
maxmin$minTGV <- as.numeric(maxmin$minTGV)
maxmin$maxTGV <- as.numeric(maxmin$maxTGV)

plotList <- list()
for (rep in 0:10)) {
  plotList[[(rep+1)]] <- ggplot(data = TGVsAll[TGVsAll$Rep==rep,], aes(x=SDGenicSt, y=zMeanGenic, group=Group, colour=scenario, linetype=scenario)) + geom_line(aes(linetype=TGVsAll$scenario[TGVsAll$Rep==rep]), size=0.5, alpha=0.4) + scale_x_reverse() +
  #geom_smooth( se=FALSE, formula=y~x+1, method="lm") + 
  xlab("Generation") + ylab("True genetic value")  + 
  scale_linetype_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        labels= c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D"), "Scenario", 
                        values=c("solid", "dotted","dashed", "dotdash", "twodash")) + 
  scale_colour_manual(breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                      labels= c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D"), "Scenario", 
                      values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1")) + 
  xlab("Genic standard deviation") + ylab("Average True Genetic Value") + guides(linetype = guide_legend(override.aes = list(size=10))) +
    geom_segment(data=maxmin[maxmin$rep==rep,], mapping=aes(x=maxGenicSD, xend=minGenicSD,
                                                            y=minTGV,  yend=maxTGV,                                    
                                                            color=scenario, linetype=scenario, group=scenario), arrow=arrow(), show.legend=FALSE, size=1.5)
}

multiplot(plotList[[1]], plotList[[2]], plotList[[3]], plotList[[7]], plotList[[8]],  cols=2)

plotListGenGain <- list()
#polot genetic gain pot replikah
for (rep in 0:10) {
  plotListGenGain[[rep+1]] <- ggplot() + geom_line(data = TGVsAll[TGVsAll$Rep==rep,], aes(x=Generation, y=zMean, group=Group, colour=scenario, linetype=scenario), size=1) + # geom_line(aes(linetype=scenario), size=0.5, alpha=0.4) + 
    xlab("Generation") + ylab("True genetic value")  + 
    scale_linetype_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), values=c("solid", "dashed", "dotted", "dotdash", "twodash"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
    scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
    xlab("Generation") + ylab("Average True Genetic Value") +
    #geom_line(data = MeanAverage, aes(x=Generation, y=MeanTGV, colour=scenario, linetype=scenario), size=1.2) + 
    ggtitle(rep) + 
    guides(group=guide_legend(nrow=6), fill=guide_legend(nrow=6), colour=guide_legend(nrow=6), linetype=guide_legend(nrow=6)) 
}

multiplot(plotListGenGain[[1]], plotListGenGain[[2]], plotListGenGain[[3]], plotListGenGain[[4]], plotListGenGain[[5]],plotListGenGain[[6]], plotListGenGain[[7]], plotListGenGain[[8]], 
          plotListGenGain[[9]], plotListGenGain[[10]],  cols=2)


#plot genska varianca po replikah
#To je plot GENSKE variance po generaijch po scenariih
plotListVar <- list()
for (rep in 0:10) {
    plotListVar[[rep+1]] <- ggplot(data = TGVsAll[TGVsAll$Rep==rep,], aes(x=Generation, y=zSdGenic, group=Group, colour=scenario, linetype=scenario)) +  geom_line(aes(linetype=scenario), size=1) + 
    xlab("Generation") + ylab("True genetic value")  + 
    scale_linetype_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                          values=c("solid", "dashed", "dotted", "dotdash", "twodash"), 
                          labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
    scale_colour_manual("Scenario", breaks = c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen"), 
                        values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), 
                        labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
    xlab("Generation") + ylab("Additive Genic Variance") + ggtitle(rep)
}

multiplot(plotListVar[[1]], plotListVar[[2]], plotListVar[[3]], plotListVar[[4]], plotListVar[[5]],plotListVar[[6]], plotListVar[[7]], plotListVar[[8]], 
          plotListVar[[9]], plotListVar[[10]],  cols=2)
