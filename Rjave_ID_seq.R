library(ROracle)

#genetski trendi pri rjavi pasmi na kodo lastnosti 164 - mleko kg

drv <- dbDriver("Oracle")

# Create the connection string
host <- "172.16.1.32"
port <- 1521
sid <- "govedo"
connect.string <- paste("(DESCRIPTION=","(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))", "(CONNECT_DATA=(SID=", sid, ")))", sep = "")


con <- dbConnect(drv, username="janao", password="job24kv5", dbname=connect.string)

poizvedba <- paste("SELECT ziv.DRZ_ORIG_ZIVAL
  || ziv.STEV_ORIG_ZIVAL ID_ZIVALI,
  ziv.ZIV_ID_SEQ,
  ziv.DAT_ROJSTVO,
  ziv.SIF_SPOL
FROM zivali ziv
WHERE ziv.SP1_SIFRA_PASMA =1")

tabela <- fetch(dbSendQuery(con,poizvedba))
tabela$DAT_ROJSTVO <- as.Date(tabela$DAT_ROJSTVO, format="%Y-%m-%d")
tabela$DAT_ROJSTVO <- format(tabela$DAT_ROJSTVO, format="%d.%m.%Y")
write.csv(tabela, "~/Genotipi/Genotipi_CODES/Rjave_seq_ID.csv", quote=F, row.names=F)
