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
and    ZGPT_ZGPVSN_SIFRA in (20,21,22,23)  
and ZGPT_VZOREC_NAROCEN_ZIV=10
ORDER BY rejec")
#tabela na Govedu, kjer so vnešeni IDji prejetih živali (Andreja)
tabela <- fetch(dbSendQuery(con,poizvedba))
tabela$Rejec <- paste(tabela$IME_REJCA, tabela$PRIIMEK_REJCA, sep=" ")

sumTabGov <- as.data.frame(table(tabela$Rejec))
colnames(sumTabGov) <- c('Rejec', 'PrejetoStVzorcevGov')

#tabela števil, prepisanih iz seznamov kontrolorjev - preverjeno s številom vzorcev
mojaTab <- read.csv('/home/jana/Documents/F4F/OdbiraZivali/StPrejetihVzorcev_Skupno.csv', na.strings='')

#preveri, ali se ujema število vzorcev po rejcu - Andrejina tabela na Govedu in moja - prepisane številke iz seznamov kontrolorjev (in preverjeno št. vzorcev)
nrow(together)
sum(together$PrejetoStVzorcev)
sum(together$PrejetoStVzorcevGov)
together[together$PrejetoStVzorcev != together$PrejetoStVzorcevGov,]

#preveri, ali so vse živali iz prvotnega seznama
zivali <- read.csv('/home/jana/Documents/F4F/OdbiraZivali/CelSeznamAplusB_15032017.csv')
nrow(tabela) == length(tabela$ZIV_ID_SEQ %in% zivali$ZIV_ID_SEQ)





##############################################################################
#Preveri ID-je vrnjenih vzorcev iz Weatherbys
###########################################################################
W_IDs <- read.csv('~/Genotipi/Genotipi_DATA/Rjava_TEMP/F4F_11072017/we_bl_sample_map_10072017.txt', sep="\t")

poizvedba <- paste("SELECT  ZGPT_ZGPVSN_SIFRA,ZGPT_ID_SEQ seq,
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
and    ZGPT_ZGPVSN_SIFRA in (20,21,22,23)
and ZGPT_VZOREC_NAROCEN_ZIV=10
ORDER BY priimek_rejca,  DRZ_ORIG_ZIVAL
|| STEV_ORIG_ZIVAL")

poizvedbaA <- paste("select * from andrejabo.anzelak")


G_IDs <- fetch(dbSendQuery(con,poizvedba))
G_IDs_A <- fetch(dbSendQuery(con,poizvedbaA)) #Anzelakovi
G_IDs_A$ZIVLJENSKA_ST <- gsub(" ","", G_IDs_A$ZIVLJENSKA_ST)

G_IDs <- c(G_IDs$ZIVAL, G_IDs_A$ZIVLJENSKA_ST) #Skupno naši poslani

length(intersect(G_IDs, W_IDs$ID)) #Koliko jih štima z Weatherbys - VSI!



###############3
#Naredi pedigre
##################

pedigreS <- paste("SELECT

  ziv.ZIV_ID_SEQ,
  ziv.ZIV_MATI_SEQ mati,
  ziv.ZIV_OCE_SEQ oce,
    ziv.DRZ_ORIG_ZIVAL
  || ziv.STEV_ORIG_ZIVAL zival,
  ziv.DRZ_STEV_ORIG_MATI || ziv.ZIV_STEV_ORIG_MATI matist,
  ziv.DRZ_STEV_ORIG_OCE || ziv.ZIV_STEV_ORIG_OCE ocest,
  ziv.DAT_ROJSTVO,
  ziv.SIF_SPOL
FROM
  govedo.zivali ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma
WHERE
  pasma.ZIV_ID_SEQ   =ziv.ZIV_ID_SEQ
AND (pasma.PASMA_SPADA=1 or pasma.pasma_spada=9)")

RjPed <- fetch(dbSendQuery(con,pedigreS))
length(intersect(RjPed$ZIVAL, W_IDs$ID)) #Koliko jih štima z Weatherbys - VSI!
RjPed$DAT_ROJSTVO <- as.Date(RjPed$DAT_ROJSTVO)
write.csv(RjPed[,c(4,5,6,7,8)], "~/Genotipi/Genotipi_DATA/Rjava_TEMP/F4F_11072017/Pedigre_11072017.csv", row.names=F, quote=F)

RjPed_F4F <- RjPed[RjPed$ZIVAL %in% G_IDs, ]
write.csv(RjPed_F4F[,c(4,5,6,7, 8)], "~/Genotipi/Genotipi_DATA/Rjava_TEMP/F4F_11072017/Pedigre_F4F.csv", row.names=F, quote=F)


