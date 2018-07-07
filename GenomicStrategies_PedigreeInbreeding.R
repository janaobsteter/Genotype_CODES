inb <- read.csv("~/InbreedingALL.csv")
colnames(inb)[c(1,2)] <- c("Generation", "F")



ScenarioDF <- data.frame()
for (strategy in unique(inb$strategy)) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen","BmGen",  "Gen")) {
    for (rep in 0:19) {
      for (interval in c(0,20,40)) {
        inb_temp <- inb[(inb$strategy==strategy) & (inb$scenario==scenario) & 
                          (inb$rep==rep) & (inb$Generation %in% interval:(interval+20)),]
        genInb <- aggregate(inb_temp$F ~ inb_temp$Generation, FUN="mean")
        colnames(genInb) <- c("Gen", "inbreeding")
        genInb$y = log(1 - genInb$inbreeding)
        genInb$t = genInb$Gen - min(genInb$Gen) + 1
        fit = MASS:::rlm(genInb$y ~ genInb$t, maxit=2000)
        # lahko bi uporabilli tudi lm(), ampak je rlm bolj robusten (pri rastlinskih simulacijah vidim precejsnje spremembe skozi cas, mislim, da pri tebi tukaj ne bi smelo biti vecjih razlik med rlm in lm)
        dF = as.data.frame(1 - exp(coef(fit)[2]))[1,1]
        Ne = 1 / (2 * dF)
        ScenarioDF$dF <- as.numeric(ScenarioDF$dF)
        ScenarioDF$Ne <- as.numeric(ScenarioDF$Ne)
        ScenarioDF$Rep <- as.numeric(ScenarioDF$Rep)
        ScenarioDF <- rbind(ScenarioDF, data.frame(Strategy=strategy, 
                                                   Scenario = scenario, Interval = interval, 
                                                   Rep = rep, dF = dF, Ne = Ne))
        colnames(ScenarioDF) <- c("Strategy", "Scenario", "Interval", "Rep", "dF", "Ne")
        
      }
    }
  }
}
