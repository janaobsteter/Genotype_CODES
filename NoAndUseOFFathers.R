fu <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/FatherUse_20032019.csv")[-1,]
fuOCS <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/FatherUse_OCS_19032019.csv")[-1,]
fuOCS1 <- read.csv("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/FatherUse_OCS_20032019_5055.csv")[-1,]
fuOCS <- rbind(fuOCS, fuOCS1)
fuOCS$Strategy <- "OCS"
table(fuOCS$Degree)
table(fu$Degree, fu$Strategy)

fu$Use <- fu$Use + 1
fuOCS$Use <- fuOCS$Use + 1


fu$Degree <- as.factor(fu$Degree)
fuOCS$Degree <- as.factor(fuOCS$Degree)
histogram(fuOCS$Use[fuOCS$Degree == 15], breaks =20)
histogram(fuOCS$Use[fuOCS$Degree == 30], breaks =20)
histogram(fuOCS$Use[fuOCS$Degree == 45], breaks =20)
histogram(fuOCS$Use[fuOCS$Degree == 50], breaks =20)
histogram(fuOCS$Use[fuOCS$Degree == 55], breaks =20)
histogram(fuOCS$Use[fuOCS$Degree == 60], breaks =20)
histogram(fuOCS$Use[fuOCS$Degree == 75], breaks =20)

fu$Degree <- revalue(fu$Degree, c("Class" = "PT", "Gen" = "GT"))

fu <- rbind(fu, fuOCS)
fu$Group <- paste0(fu$Strategy, fu$Degree)
table(fu$Group)
fu <- fu[fu$Group %in% c("SU55PT", "SU55GT", "SU51GT","OCS15","OCS30","OCS45","OCS50", "OCS55","OCS60","OCS75" ),]
fu$Group <- factor(fu$Group, level=c("SU55PT","SU55GT","SU51GT","OCS15","OCS30","OCS45","OCS50", "OCS55","OCS60","OCS75"))


fuA <- summarySE(data=fu, groupvars = c("Group", "Rep"), measurevar = "Use")
summarySE(data=fuA, groupvars = c("Group"), measurevar = "Use")

colnames(fu)[1] <- "Scenario"

fu$Group <- paste0(fu$Strategy, fu$Scenario)


model <- lm(Use ~ Group, data=fuA)
marginal = emmeans(model, ~ Group)
pairs(marginal,
      adjust="tukey")
CLD = cld(marginal, sort=FALSE,
          alpha   = 0.05,
          Letters = letters, adjust="tukey") 
CLD


#Number of fathers

NO <- read.table("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/NoFather_Trunc_20032019.txt")[-1,]
NO_OCS <- read.table("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/NoFather_OCS_11022019.txt")[-1,]
NO_OCS1 <- read.table("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/NoFather_OCS_20032019_5055.txt")[-1,]
NO_OCS <- rbind(NO_OCS, NO_OCS1)
table(NO_OCS$Degree)

NO <- NO[NO$Generation != 40,]
NO_OCS <- NO_OCS[NO_OCS$Generation != 40,]

NO_OCS$Strategy <- "OCS"
colnames(NO_OCS)[4] <- "Scenario"
NO$Scenario <- as.factor(NO$Scenario)
NO_OCS$Scenario <- as.factor(NO_OCS$Scenario)
NO$Scenario <- revalue(NO$Scenario, c("Class" = "PT", "Gen" = "GT"))

NO <- rbind(NO, NO_OCS)
NO$Group <- paste0(NO$Strategy, NO$Scenario)
table(NO$Group)
NO <- NO[NO$Group %in% c("SU55PT", "SU55GT", "SU51GT","OCS15","OCS30","OCS45","OCS50", "OCS55","OCS60","OCS75" ),]
NO$Group <- factor(NO$Group, level=c("SU55PT","SU55GT","SU51GT","OCS15","OCS30","OCS45","OCS50", "OCS55","OCS60","OCS75"))

#povprečje povprečij
NOa <- summarySE(data=NO, groupvars = c("Group", "Rep"), measurevar = "NoFathers")
summarySE(data=NOa, groupvars = c("Group"), measurevar = "NoFathers")


model <- lm(NoFathers ~ Group, data=NOa)
marginal = emmeans(model, ~ Group)
pairs(marginal,
      adjust="tukey")
CLD = cld(marginal, sort=FALSE,
          alpha   = 0.05,
          Letters = letters, adjust="tukey") 
CLD



ped <- read.csv("CriterionandTGV.csv")
ped$Sel <- ifelse(ped$nMating > 1, TRUE, FALSE)
table(ped$Sel)
table(ped$sex)
ped <- ped[ped$sex == "M",]
head(ped)
cor(ped$gvNormUnres1, ped$Criterion)
ggplot(data=ped, aes(x=gvNormUnres1, y=Criterion, colour=Sel)) + geom_point()
