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

compare <- merge(MCPklasika, MCPspek, by="ID", all.x=TRUE)
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

#ujemanje RCT na vzorcih, ki so bili v istem mesecu narejeni na klasiki in spektrofotometrično
rct <- melt(SAME, id.vars = "ID", measure.vars = c("rct_sp", "rct_kl") )
rct <- rct[rct$value != 0,]
ggplot(rct, aes(x=ID, y=value, group=variable, fill=variable)) + geom_histogram(stat="identity", position="identity") + 
  xlab("ID živali") + ylab("a30") + scale_fill_manual("Meritev", breaks=c("a30_sp", "a30_kl"), values=c("steelblue2", "darkblue"), labels=c("Klasično", "Spektrofotometrično")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_text(size=14), axis.text.y = element_text(size=14), 
        legend.position = "bottom", legend.title = element_text(size=14), legend.text = element_text(size=14))

#ujemanje zimske in pomlade analize
oboje <- intersect(compare$ID[compare$Sezona_kl=="Zima"], compare$ID[compare$Sezona_kl=="Pomlad"])
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
