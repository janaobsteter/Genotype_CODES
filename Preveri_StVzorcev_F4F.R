drv <- dbDriver("Oracle")

# Create the connection string
host <- "172.16.1.32"
port <- 1521
sid <- "govedo"
connect.string <- paste("(DESCRIPTION=","(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))", "(CONNECT_DATA=(SID=", sid, ")))", sep = "")


con <- dbConnect(drv, username="janao", password="job24kv5", dbname=connect.string)

poizvedba <- paste("SELECT  ZGPT_ZGPVSN_SIFRA,ZGPT_ID_SEQ seq, zivali.ziv_id_seq,
   ( DRZ_ORIG_ZIVAL
   || STEV_ORIG_ZIVAL) zival,
 govedo_splosne_procedure.poisci_ime (ZGPT_ZIV_ID_SEQ) ime,
 govedo_splosne_procedure. preberi_kratko_ime_pasma ( govedo_splosne_procedure.preberi_sifro_pasme (ZGPT_ZIV_ID_SEQ)) pasma,
 govedo_splosne_procedure. preberi_naslov_rejca_1vrsta (ZGPT_LOKACIJA_REJEC) rejec,
 GOVEDO_SPLOSNE_PROCEDURE.POISCI_KMG_MID (ZGPT_LOKACIJA_REJEC) kmgmid,
 GOVEDO_SPLOSNE_PROCEDURE.preberi_ime_rejca (ZGPT_LOKACIJA_REJEC) ime_rejca,
 GOVEDO_SPLOSNE_PROCEDURE.preberi_priimek_rejca (ZGPT_LOKACIJA_REJEC) priimek_rejca,
 govedo_splosne_procedure.preberi_ime_kontrolorja (ZGPT_KONTROLOR) kontrolor_ime,
 ZGPT_KONTROLOR AS kontrolor_sif
FROM govedo.ZIVALI_GP_TEST, govedo.zivali
WHERE ZGPT_ZIV_ID_SEQ = ZIV_ID_SEQ
and    ZGPT_ZGPVSN_SIFRA in (20,21)  
and ZGPT_VZOREC_NAROCEN_ZIV=10
ORDER BY rejec")
#tabela na Govedu, kjer so vnešeni IDji prejetih živali (Andreja)
tabela <- fetch(dbSendQuery(con,poizvedba))
tabela$Rejec <- paste(tabela$IME_REJCA, tabela$PRIIMEK_REJCA, sep=" ")

sumTabGov <- as.data.frame(table(tabela$Rejec))
colnames(sumTabGov) <- c('Rejec', 'PrejetoStVzorcevGov')

#tabela števil, prepisanih iz seznamov kontrolorjev - preverjeno s številom vzorcev
mojaTab <- read.csv('/home/jana/Documents/F4F/OdbiraZivali/StPrejetihVzorcev_16052017.csv', na.strings='')

#preveri, ali se ujema število vzorcev po rejcu - Andrejina tabela na Govedu in moja - prepisane številke iz seznamov kontrolorjev (in preverjeno št. vzorcev)
together <- merge(mojaTab, sumTabGov, by='Rejec', all = T)
nrow(together)
sum(together$PrejetoStVzorcev)
sum(together$PrejetoStVzorcevGov)
together[together$PrejetoStVzorcev != together$PrejetoStVzorcevGov,]

#preveri, ali so vse živali iz prvotnega seznama
zivali <- read.csv('/home/jana/Documents/F4F/OdbiraZivali/CelSeznamAplusB_15032017.csv')
nrow(tabela) == length(tabela$ZIV_ID_SEQ %in% zivali$ZIV_ID_SEQ)

