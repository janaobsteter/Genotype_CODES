library(ggplot2)
p = ggplot(data=TGVsAll, mapping=aes(x=zSdGenic, y=zMean, color=scenario, linetype=scenario)) +  geom_path(size=1/5, alpha=.5) + 
#  # Conv, PYT, Head, TwoPartTS, TwoPartTS+, TwoPartOCS  
scale_color_manual(values=c("black", "red", "green", "blue", "yellow"), name="") +  
scale_linetype_manual(values=c(      1,    1,    1,    1,    2,    1), name="") +  
xlab("Genic standard deviation") + ylab("Genetic mean") +  
scale_x_reverse(sec.axis=sec_axis(trans=~1-.,                                   
                                    name="Converted/Lost genic standard deviation")) +  
theme_bw() +  theme(legend.position="top") +  
guides(colour=guide_legend(nrow=1)) +  
geom_segment(data=tmp5, mapping=aes(x=zSdGenicMax, xend=zSdGenicMin,                                      
                                     y=zMeanStartSdGenic, yend=zMeanEndSdGenic,                                      
                                     color=scenario3, linetype=scenario3)) +  
geom_segment(data=tmp5, mapping=aes(x=zSdGenicMax, xend=zSdGenicMin,                                      
                                      y=zMeanStartSdGenic, yend=zMeanEndSdGenic,                                      
                                      color=scenario3, linetype=scenario3),               
               arrow=arrow(), show.legend=FALSE)



efficiency = function(x, f) {
  ret = NA
  if (nrow(x) > 5) {
    x = arrange_(x, "Generation")
    # fit = lm(formula=f, data=x)
    fit = MASS:::rlm(formula=f, data=x, maxit=1000)
    # fit = MASS:::lqs(formula=f, data=x)
    ret = -coef(fit)[2]
  }
  ret
}
datByRunStage = dat %>%
  group_by(run, rep, scenario, sel, self, size, scaled, crit, nCycles, target, stage) %>%
  nest()
data <- TGVsAll
datByRunStage = datByRunStage %>%
  mutate(efficiencyG=map(data, efficiency, f=zMean      ~ zSdGenic),
         efficiencyGenic=map(data, efficiency, f=zMeanGenic ~ zSdGenic))
# zMean je standardiziran genetski napredek (TBV - mean(TBV_gen_start)) / sd(TBV_gen_start)
# zMeanGenic je standardiziran genetski napredek (TBV - mean(TBV_gen_start)) / sqrt(2*sum(p*q*alpha)_gen_start


################
#Funkcija summarySE
#####################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
######################################################################

TGVsAll <- data.frame()


for (scenario in c("Class", "GenSLO", "GenSplosnaPop", "GenSLO_BmGen", "Gen")) {
  TGVs <- data.frame(Generation=40:60)
  dir = paste0('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_', scenario, '/SimulatedData/')
  ped <- read.table(paste0(dir,'/PedigreeAndGeneticValues_cat.txt'), header=T)
  ped <- ped[ped$Generation %in% 40:60,]
  TGV <- summarySE(ped, measurevar = "gvNormUnres1", groupvars = "Generation")[,c(1,3,4)]
  colnames(TGV)[1] <- c("Generation")
  TGV$zMean <- (TGV$gvNormUnres1 - TGV$gvNormUnres1[1]) / TGV$sd[1]
  TGVs <- merge(TGVs, TGV, by="Generation")
  Var <- read.table(paste0(dir,'/TotalGenicAndGeneticVariancesPerGeneration.txt'), header=T)
  Var <- Var[Var$QtnModel==1,c(1,3)]
  TGVs <- merge(TGVs, Var, by="Generation")
  TGVs$zMeanGenic <- (TGVs$gvNormUnres1 - TGVs$gvNormUnres1[1]) / TGVs$AdditGenicVar1[1]
  #TGVsAll$zSdGenic <- (sqrt(TGVsAll$AdditGenicVar1) - sqrt(TGVsAll$))
  TGVs$scenario <- scenario
  #colnames(TGVs) < c("Generation", paste0("TGV_mean", scenario), paste0("TGV_sd", scenario), paste0("zMean_", scenario), paste0("GenicVar_", scenario), paste0("zMeanGenic_", scenario))
  TGVsAll <- rbind(TGVsAll, TGVs)
}

#TGVsAll$AdditGenicVar1 <- as.numeric(TGVsAll$AdditGenicVar1)
TGVsAll$zSdGenic <- (sqrt(TGVsAll$AdditGenicVar1)) 

library(MASS)
library(ggplot2)
lm <- ggplot(data = TGVsAll, aes(x=TGVsAll$zSdGenic, y=TGVsAll$zMeanGenic, colour=TGVsAll$scenario)) + geom_path() + scale_x_reverse() + 
  geom_smooth(method='lm', se=FALSE) + 
  scale_color_hue("Shema", labels=c("Conventional", "GenomicSLO", "GenBulls on Other Cows", "GenBulls on Bull Dams", "GenBulls on All Cows")) + 
  xlab("Genic sd") + ylab("Mean genetic gain") + ggtitle("Genic variance standardisation")
lm2 <- ggplot(data = TGVsAll, aes(x=TGVsAll$zSdGenic, y=TGVsAll$zMean, colour=TGVsAll$scenario)) + geom_path() + scale_x_reverse() + 
  geom_smooth(method='lm', se=FALSE) + 
  scale_color_hue("Shema", labels=c("Conventional", "GenomicSLO", "GenBulls on Other Cows", "GenBulls on Bull Dams", "GenBulls on All Cows")) + 
  xlab("Genic sd") + ylab("Mean genetic gain") + ggtitle("Genetic variance standardisation")
rlm <- ggplot(data = TGVsAll, aes(x=TGVsAll$zSdGenic, y=TGVsAll$zMeanGenic, colour=TGVsAll$scenario)) + geom_path() + scale_x_reverse() + geom_smooth(method='rlm', se=FALSE)
library(Rmisc)
multiplot(lm, rlm, cols=2)
multiplot(lm, lm2, cols=2)
