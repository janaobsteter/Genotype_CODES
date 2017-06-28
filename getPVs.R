library(ROracle)

#genetski trendi pri rjavi pasmi na kodo lastnosti 164 - mleko kg

drv <- dbDriver("Oracle")

# Create the connection string
host <- "172.16.1.32"
port <- 1521
sid <- "govedo"
connect.string <- paste("(DESCRIPTION=","(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))", "(CONNECT_DATA=(SID=", sid, ")))", sep = "")


con <- dbConnect(drv, username="janao", password="job24kv5", dbname=connect.string)

poizvedba <- paste("SELECT ziv.ZIV_ID_SEQ, ziv.SIF_SPOL,pv.VREDNOST_12_PV, extract(year from pv.DAT_OCENA_PV)
FROM govedo.zivali ziv,
ARHIV.PLEMENSKE_VREDNOSTI pv
WHERE ziv.SP1_SIFRA_PASMA = 1
AND ziv.ZIV_ID_SEQ        =pv.PV_ZIV_ID_SEQ
AND pv.SIFRA_LAST         =164
AND ziv.DRZ_ORIG_ZIVAL    ='SI'")

tabela <- fetch(dbSendQuery(con,poizvedba))
means <- aggregate(tabela$VREDNOST_12_PV ~ tabela$`EXTRACT(YEARFROMPV.DAT_OCENA_PV)`, FUN=mean)
