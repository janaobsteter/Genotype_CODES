mcp <- read.csv("~/Documents/F4F/Rezultati_MCpSpektrofotimetrija.csv")
mcp$Date <- as.Date(mcp$DAT_ANA, format="%d-%m-%y")

mcp$RCT <- as.numeric(mcp$RCT)
mcp$A30 <- as.numeric(mcp$A30)
mcp$RCTs <- mcp$RCT * 60 
par(mfrow=c(2,1))
hist(mcp$RCTs, breaks=50)
hist(MCPA$r)
hist(mcp$A30, breaks=50)
hist(MCPA$a30.mm.)

library(ggplot2)
library(plyr)
mcp$SP1_SIFRA_PASMA <- revalue(factor(mcp$SP1_SIFRA_PASMA), c("1" = "Rjava", "2" = "Lisasta", "3" = "Črno bela", "5" = "Cika"))
ggplot(data=mcp[mcp$SP1_SIFRA_PASMA %in% c("Rjava", "Lisasta", "Črno bela"),], aes(x = SP1_SIFRA_PASMA, y = RCT, fill=SP1_SIFRA_PASMA)) + stat_summary(fun.y="mean", geom="bar") + 
  scale_fill_discrete("PASMA") + xlab("PASMA") + theme(legend.position = "none", axis.title =  element_text(size=14), axis.text = element_text(size=14), legend.text = element_text(size=12)) 
ggplot(data=mcp[mcp$SP1_SIFRA_PASMA %in% c("Rjava", "Lisasta", "Črno bela"),], aes(x = SP1_SIFRA_PASMA, y = A30, fill=SP1_SIFRA_PASMA)) + stat_summary(fun.y="mean", geom="bar") + 
  scale_fill_discrete("PASMA") + xlab("PASMA")  + theme(legend.position = "none", axis.title =  element_text(size=14), axis.text = element_text(size=14), legend.text = element_text(size=12)) 
scale_fill_discrete("PASMA") + xlab("PASMA") 
ggplot(data=mcp[mcp$SP1_SIFRA_PASMA %in% c("Rjava", "Lisasta", "Črno bela") & mcp$RCT != 0 & mcp$BELJ > 2 & mcp$BELJ < 5.5,], aes(x = BELJ, y = RCT)) + geom_point(alpha=0.1) + geom_smooth()  + ylim(c(0, 10))
ggplot(data=mcp[mcp$SP1_SIFRA_PASMA %in% c("Rjava", "Lisasta", "Črno bela") & mcp$RCT != 0 & mcp$MASC > 2 & mcp$MASC < 5.5,], aes(x = MASC, y = RCT)) + geom_point() + geom_smooth() 
ggplot(data=mcp[mcp$SP1_SIFRA_PASMA %in% c("Rjava", "Lisasta", "Črno bela") & mcp$A30 != 0 & mcp$BELJ > 2 & mcp$BELJ < 5.5,], aes(x = BELJ, y = A30)) + geom_point() + geom_smooth() 

  scale_fill_discrete("PASMA")

aggregate(mcp$A30 ~ mcp$SP1_SIFRA_PASMA, FUN="mean")
mcp[mcp$SP1_SIFRA_PASMA == 0,]

summarySE(data=mcp, measurevar = "A30", groupvars = "SP1_SIFRA_PASMA")

mcp$SP1_SIFRA_PASMA <- as.factor(mcp$SP1_SIFRA_PASMA)
library(ggplot2)
ggplot(data=mcp[mcp$SP1_SIFRA_PASMA %in% 1:3,], aes(x=A30, fill=SP1_SIFRA_PASMA)) + geom_density()

#preveri skladnost klasike in spektrofotometrije
`%not_in%` <- purrr::negate(`%in%`)
length(intersect(MCPA$ID, mcp$ID_ZIVALI))
sum((unique(MCPA$ID) %not_in% unique(mcp$ID_ZIVALI)))
unique(MCPA$ID)[((unique(MCPA$ID) %not_in% unique(mcp$ID_ZIVALI)))]
sum((unique(MCPA$ID) %in% unique(mcp$ID_ZIVALI)))
length(unique(MCPA$ID))

MCPklasika <- MCPA[,c("ID", "r.s.", "k20.s.", "a30.mm.", "Mascoba", "proteini", "Season", "Datum")]
colnames(MCPklasika) <- c("ID", "rct_kl", "k20_kl", "a30_kl", "masc_kl", "belj_kl", "Sezona_kl", "Datum_kl")       
MCPspek <- mcp[,c("ID_ZIVALI", "SP1_SIFRA_PASMA", "Date", "RCT", "A30", "BELJ", "MASC")]
colnames(MCPspek) <- c("ID", "PASMA", "Datum_sp", "rct_sp",  "a30_sp", "belj_sp", "masc_sp")       

MCPklasika$Year <- format(MCPklasika$Datum_kl, "%y")
MCPspek$Year <- format(MCPspek$Datum_sp, "%y")
compare <- merge(MCPklasika, MCPspek, by=c("ID", "Year"), all.x=TRUE)
cor(compare$rct_kl[compare$Year == 17], compare$rct_sp[compare$Year == 17], use="pairwise.complete.obs")
cor(compare$rct_kl[compare$Year == 18], compare$rct_sp[compare$Year == 18], use="pairwise.complete.obs")
cor(compare$a30_kl[compare$Year == 17], compare$a30_sp[compare$Year == 17], use="pairwise.complete.obs")
cor(compare$a30_kl[compare$Year == 18], compare$a30_sp[compare$Year == 18], use="pairwise.complete.obs")
aggregate(compare$a30_kl ~ compare$Year, FUN="mean")
aggregate(compare$a30_sp ~ compare$Year, FUN="mean")
aggregate(compare$rct_kl ~ compare$Year, FUN="mean")
aggregate(compare$rct_sp ~ compare$Year, FUN="mean")

compare <- merge(MCPklasika, MCPspek, by=c("ID"), all.x=TRUE)

compare <- compare[,c("ID", "PASMA", "Sezona_kl", "Datum_kl", "Datum_sp", "rct_kl", "rct_sp", "a30_kl", "a30_sp", "k20_kl", "masc_kl", "masc_sp", "belj_kl", "belj_sp")]
compare$MonthKl <- format(compare$Datum_kl, "%m-%y")
compare$MonthSp <- format(compare$Datum_sp, "%m-%y")
compare <- compare[order(compare$ID, compare$Datum_kl, compare$Datum_sp),]

cor(compare$rct_kl, compare$rct_sp, use="pairwise")
cor(compare$a30_kl, compare$a30_sp, use="pairwise")

#ujemanje a30 na vzorcih, ki so bili v istem mesecu narejeni na klasiki in spektrofotometrično
SAME <- compare[compare$MonthKl == compare$MonthSp,]

library(reshape)
a30 <- melt(SAME, id.vars = "ID", measure.vars = c("a30_sp", "a30_kl") )
a30 <- a30[a30$value != 0,]
ggplot(a30, aes(x=ID, y=value, group=variable, fill=variable)) + geom_histogram(stat="identity", position="identity") + 
  xlab("ID živali") + ylab("a30") + scale_fill_manual("Meritev", breaks=c("a30_kl", "a30_sp"), values=c("steelblue2", "darkblue"), labels=c("Klasično", "Spektrofotometrično")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_text(size=14), axis.text.y = element_text(size=14), 
        legend.position = "bottom", legend.title = element_text(size=14), legend.text = element_text(size=14))
cor(SAME$a30_kl, SAME$a30_sp, use="pairwise.complete.obs")
model <- lm(SAME$a30_kl ~ SAME$a30_sp)
mean(model$residuals^2)
mean(SAME$a30_kl, na.rm=TRUE)
mean(SAME$a30_sp, na.rm = TRUE)

#ujemanje RCT na vzorcih, ki so bili v istem mesecu narejeni na klasiki in spektrofotometrično
rct <- melt(SAME, id.vars = "ID", measure.vars = c("rct_sp", "rct_kl") )
rct <- rct[rct$value != 0,]
ggplot(rct, aes(x=ID, y=value, group=variable, fill=variable)) + geom_histogram(stat="identity", position="identity") + 
  xlab("ID živali") + ylab("a30") + scale_fill_manual("Meritev", breaks=c("a30_sp", "a30_kl"), values=c("steelblue2", "darkblue"), labels=c("Klasično", "Spektrofotometrično")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_text(size=14), axis.text.y = element_text(size=14), 
        legend.position = "bottom", legend.title = element_text(size=14), legend.text = element_text(size=14))

#ujemanje zimske in pomlade analize
#klasika
a30_kl_cor <- merge(MCPklasika[MCPklasika$Year == 17, c("ID", "a30_kl")], MCPklasika[MCPklasika$Year == 18, c("ID", "a30_kl")], by="ID")
cor(a30_kl_cor$a30_kl.x, a30_kl_cor$a30_kl.y, use = "pairwise.complete.obs")
rct_kl_cor <- merge(MCPklasika[MCPklasika$Year == 17, c("ID", "rct_kl")], MCPklasika[MCPklasika$Year == 18, c("ID", "rct_kl")], by="ID")
cor(rct_kl_cor$rct_kl.x, rct_kl_cor$rct_kl.y, use = "pairwise.complete.obs")

#spektrofotometrija
a30_sp_cor <- merge(MCPspek[MCPspek$Year == 17, c("ID", "a30_sp")], MCPspek[MCPspek$Year == 18, c("ID", "a30_sp")], by="ID")
cor(a30_sp_cor$a30_sp.x, a30_sp_cor$a30_sp.y, use = "pairwise.complete.obs")
rct_sp_cor <- merge(MCPspek[MCPspek$Year == 17, c("ID", "rct_sp")], MCPspek[MCPspek$Year == 18, c("ID", "rct_sp")], by="ID")
cor(rct_sp_cor$rct_sp.x, rct_sp_cor$rct_sp.y, use = "pairwise.complete.obs")

oboje <- intersect(compare$ID[compare$Sezona_kl=="Zima"], compare$ID[compare$Sezona_kl=="Pomlad"])
oboje <- intersect(compare$ID[compare$Year==17], compare$ID[compare$Year==18])
sezona  <- unique(melt(compare[compare$ID %in% oboje & compare$a30_kl != 0,], id.vars = c("ID", "Sezona_kl"), measure.vars = c("a30_kl") ))
library(tidyr)
sezonaS <- spread(data=sezona,  key = Sezona_kl, value=value)
cor(sezonaS$Pomlad, sezonaS$Zima, use = "pairwise.complete.obs")
#rank correlation
cor.test( ~ Zima + Pomlad, 
          data=sezonaS,
          method = "spearman",
          continuity = FALSE,
          conf.level = 0.95)
ggplot(sezona, aes(x=ID, y=value, group=Sezona_kl, fill=Sezona_kl, colour=Sezona_kl)) + geom_path() + #geom_histogram(stat="identity", position="identity") + 
  xlab("ID živali") + ylab("a30") + scale_fill_manual("Meritev", breaks=c("a30_sp", "a30_kl"), values=c("steelblue2", "darkblue"), labels=c("Klasično", "Spektrofotometrično")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_text(size=14), axis.text.y = element_text(size=14), 
        legend.position = "bottom", legend.title = element_text(size=14), legend.text = element_text(size=14))

ggplot(sezonaS, aes(x=Zima, y=Pomlad)) + geom_point() + #geom_histogram(stat="identity", position="identity") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_text(size=14), axis.text.y = element_text(size=14), 
        legend.position = "bottom", legend.title = element_text(size=14), legend.text = element_text(size=14))

