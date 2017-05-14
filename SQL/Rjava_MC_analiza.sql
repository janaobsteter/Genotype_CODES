select distinct mc.RM_SIF_LOKACIJA creda from GOVEDO.REJCI_MLEKARNA mc, govedo.zivali ziv where ziv.CRE_SIFRA_CREDA=mc.RM_SIF_LOKACIJA and mc.RM_REJEC_ODDAJA_MLEKO=1;




--rejci, ki oddajajo mleko in so v AT kontroli
SELECT DISTINCT
  mc.RM_SIF_LOKACIJA creda,
  lok.*
FROM
  GOVEDO.REJCI_MLEKARNA mc,
  GOVEDO.LOKACIJE_STEVILO_ZIVALI lok,
  GOVEDO.STATUSI_ZIVALI_CREDA stat
  WHERE
  lok.LOKACIJA              =mc.RM_SIF_LOKACIJA
  and mc.RM_SIF_LOKACIJA=stat.CRE_SIFRA_CREDA
  and stat.SIF_STATUS_CREDA=5
  AND mc.RM_REJEC_ODDAJA_MLEKO=1
  and lok.DATUM_POROCILA='31.01.2017';
  
  
--rejci, ki oddajajo mleko in so v AT kontroli in imajo rjave krave
SELECT DISTINCT
  mc.RM_SIF_LOKACIJA creda,
  lok.*
FROM
  GOVEDO.REJCI_MLEKARNA mc,
  GOVEDO.LOKACIJE_STEVILO_ZIVALI lok,
  GOVEDO.STATUSI_ZIVALI_CREDA stat
  WHERE
  lok.LOKACIJA              =mc.RM_SIF_LOKACIJA
  and mc.RM_SIF_LOKACIJA=stat.CRE_SIFRA_CREDA
  and stat.SIF_STATUS_CREDA=5
  AND mc.RM_REJEC_ODDAJA_MLEKO=1
  and lok.DATUM_POROCILA='30.04.2017'
  and lok.RJAR not in 0
;
  

--rejci, ki oddajajo mleko in so v AT kontroli in število rjavih živali v laktaciji --> SSI rjavih živali

drop table Reje_MC_11052017 ;
create table Reje_MC_11052017 as;
SELECT DISTINCT
  zivali_mc.creda,
  lok.RJAR RjKrave,
  round(avg(pv.VREDNOST_12_PV),2) avgSSI,
  round(stddev(pv.VREDNOST_12_PV),3) SSIsd,
  count(distinct zivali_mc.ziv_id_seq) stZivalisPV
FROM
  ARHIV.PLEMENSKE_VREDNOSTI pv,
  --select aktivne zivali na podrocju mlekarne, ki oddajajo mleko in so v AT kontroli
  (SELECT DISTINCT
  mc.RM_SIF_LOKACIJA creda,
  ziv.ziv_id_seq
FROM
  GOVEDO.REJCI_MLEKARNA mc,
  GOVEDO.STATUSI_ZIVALI_CREDA stat,
  govedo.zivali ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma
WHERE
  ziv.CRE_SIFRA_CREDA       =mc.RM_SIF_LOKACIJA
AND mc.RM_SIF_LOKACIJA      =stat.CRE_SIFRA_CREDA
AND ziv.ZIV_ID_SEQ          =pasma.ZIV_ID_SEQ
AND stat.SIF_STATUS_CREDA   =5
AND mc.RM_REJEC_ODDAJA_MLEKO=1
AND ziv.AKTIVNA             =1
and pasma.PASMA_SPADA=1
and ziv.SIF_SPOL=2
) zivali_mc,
GOVEDO.LOKACIJE_STEVILO_ZIVALI lok
  WHERE
  lok.LOKACIJA=zivali_mc.creda
and zivali_mc.ziv_id_seq=pv.PV_ZIV_ID_SEQ(+)
  and (extract(year from pv.DAT_OCENA_PV)=2016)
  and (extract(month from pv.DAT_OCENA_PV)=11) 
 and pv.SIFRA_LAST(+)=79
  and lok.RJAR not in 0
  --and lok.RJar >= 10 --več kot 10 rjavih krav
  and lok.DATUM_POROCILA='30.04.2017'
  group by zivali_mc.creda, lok.RJAR;
  
  
  
--rejci, ki oddajajo mleko in so v AT kontroli (in imajo več kot 10 rajvih krav) --> živali - ali so v laktaciji ali ne
SELECT ziv.ZIV_ID_SEQ, ziv.IME_ZIVAL,ziv.SIF_BREJA ,ziv.CRE_SIFRA_CREDA,ziv.DRZ_ORIG_ZIVAL || ziv.STEV_ORIG_ZIVAL ID_zivali, 
ziv.DAT_ROJSTVO, ziv.ZIV_MATI_SEQ mati, ziv.ZIV_OCE_SEQ oce, ziv.DRZ_STEV_ORIG_MATI || ziv.ZIV_STEV_ORIG_MATI matiID, ziv.DRZ_STEV_ORIG_OCE || ziv.ZIV_STEV_ORIG_OCE oceID
,max(tel.ZAP_TELITEV) StTel, max(tel.DAT_TELITEV) ZadnjaTel, pv.VREDNOST_12_PV
FROM
Reje_MC_11052017 reje, GOVEDO.ZIVALI ziv,
GOVEDO.ZIVALI_PASMA_SPADA pasma,
GOVEDO.TELITVE tel, GOVEDO.LAKTACIJE lak,
ARHIV.PLEMENSKE_VREDNOSTI pv
  WHERE
  ziv.ZIV_ID_SEQ=pv.PV_ZIV_ID_SEQ(+)
  and pv.SIFRA_LAST(+)=79
  and ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ 
  and tel.TEL_ID_SEQ=lak.LAK_TEL_ID_SEQ and
  ziv.CRE_SIFRA_CREDA=reje.CREDA
  and pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
  and ziv.SIF_SPOL=2
  and ziv.AKTIVNA=1
  and pasma.PASMA_SPADA=1
 -- and ziv.CRE_SIFRA_CREDA=32343
 --group by ziv.CRE_SIFRA_CREDA
--group by ziv.ZIV_ID_SEQ
group by ziv.ZIV_ID_SEQ, ziv.DRZ_ORIG_ZIVAL,ziv.STEV_ORIG_ZIVAL, ziv.DAT_ROJSTVO, ziv.ZIV_MATI_SEQ, ziv.ZIV_OCE_SEQ,
 ziv.DRZ_STEV_ORIG_MATI, ziv.ZIV_STEV_ORIG_MATI , ziv.DRZ_STEV_ORIG_OCE , ziv.ZIV_STEV_ORIG_OCE, ziv.sif_breja, ziv.CRE_SIFRA_CREDA, pv.VREDNOST_12_PV, ziv.IME_ZIVAL
;

select * from GOVEDO.LOKACIJE_STEVILO_ZIVALI lok where lok.LOKACIJA=32343 and lok.DATUM_POROCILA='31.01.2017';
  
select * from ARHIV.PLEMENSKE_VREDNOSTI pv where pv.PV_ZIV_ID_SEQ=2721347 and pv.SIFRA_LAST=79;
  
  
--najbolše reje glede na SSI mleko
SELECT DISTINCT
  govedo_splosne_procedure.PREBERI_NASLOV_REJCA_3VRSTE (ziv.CRE_SIFRA_CREDA),
  ziv.CRE_SIFRA_CREDA
FROM
  zivali ziv
WHERE
  ziv.CRE_SIFRA_CREDA IN (8954, 3295, 9072);
  
--najslabše reje glede na SSI
  SELECT DISTINCT
  govedo_splosne_procedure.PREBERI_NASLOV_REJCA_3VRSTE (lok.LOKACIJA),
  lok.RJAR RjKrave
FROM
  GOVEDO.LOKACIJE_STEVILO_ZIVALI lok
WHERE
  lok.LOKACIJA IN (1868, 1626, 31992, 2289, 10176 ,8954, 482, 15731, 3295, 32162,
  2091, 9072 ,4219, 7150, 8619 ) and lok.DATUM_POROCILA='31.01.2017';
  
  
  
--črede z več kot 20 rjavimi kravami
  SELECT DISTINCT
  govedo_splosne_procedure.PREBERI_NASLOV_REJCA_3VRSTE (lok.LOKACIJA),
  lok.RJAR RjKrave, lok.lokacija
FROM
  GOVEDO.LOKACIJE_STEVILO_ZIVALI lok
WHERE
  lok.LOKACIJA IN (9072,  8180,   150, 32717 , 4128,  6658, 13035,  8060 ,  143 , 7480 , 4950  , 178  ,3989,  1171 ,32749 ,15631 , 9570 ) and lok.DATUM_POROCILA='31.01.2017';

--poglej po očetih - koliko aktivnih hčera imajo na območju MC
SELECT
  COUNT (DISTINCT ziv.ZIV_ID_SEQ),
  ziv.ZIV_OCE_SEQ seqOce,
  gen.GEN_CHIP,
  pasma.PASMA_SPADA
FROM
  govedo.zivali ziv,
  GOVEDO.REJCI_MLEKARNA rm,
  GENOTIPIZIRANE_ZIVALI gen,
  GOVEDO.ZIVALI_PASMA_SPADA pasma
WHERE
  RM_SIF_LOKACIJA         = CRE_SIFRA_CREDA
AND gen.ZIV_ID_SEQ(+)     =ziv.ZIV_OCE_SEQ
AND pasma.ZIV_ID_SEQ      =ziv.ZIV_OCE_SEQ
AND RM_REJEC_ODDAJA_MLEKO = 1
AND ziv.AKTIVNA           =1
AND ziv.SIF_SPOL          =2
AND pasma.PASMA_SPADA     =1
GROUP BY
  ziv.ZIV_OCE_SEQ,
  gen.GEN_CHIP,
  pasma.PASMA_SPADA;


--očetje aktivnih ženskih živali na območju MC
SELECT
  COUNT(DISTINCT mc_zivali.seqZiv) StPotomk,
  OceZiv.ZIV_ID_SEQ OceSeq,
  gen.GEN_CHIP,
  OceZiv.DAT_ROJSTVO OceDat,
  genRj.GEN_KAPPA_KAZEIN kappaGen,
  genRj.GEN_BETA_GENOTIP betaGen,
  gp.ZGP_VREDNOST
FROM
  zivali OceZiv,
  (
    SELECT
      ziv.ZIV_ID_SEQ seqZiv,
      ziv.ZIV_OCE_SEQ seqOce,
      ziv.DAT_ROJSTVO
    FROM
      govedo.zivali ziv,
      GOVEDO.REJCI_MLEKARNA rm
    WHERE
      RM_SIF_LOKACIJA         = CRE_SIFRA_CREDA
    AND RM_REJEC_ODDAJA_MLEKO = 1
    AND ziv.aktivna           =1
    AND ziv.sif_spol          =2
  )
  mc_zivali,
  GOVEDO.ZIVALI_PASMA_SPADA pasma,
  GENOTIPIZIRANE_ZIVALI gen,
  govedo.zivali_gp gp,
  GOVEDO.GEN_RJAVA genRj
WHERE
genRj.GEN_ZIV_ID_SEQ(+)=OceZiv.ZIV_ID_SEQ
and gp.ZGP_ZIV_ID_SEQ(+)=OceZiv.ZIV_ID_SEQ and
gp.ZGP_SGPL_SIFRA(+) = 100
and  mc_zivali.seqOce   =OceZiv.ziv_id_seq
  and gen.ZIV_ID_SEQ(+)=OceZiv.ZIV_ID_SEQ
AND pasma.ZIV_ID_SEQ =OceZiv.ZIV_ID_SEQ
AND pasma.PASMA_SPADA=1 group by   
OceZiv.ZIV_ID_SEQ,
  gen.GEN_CHIP,
  OceZiv.DAT_ROJSTVO ,
  genRj.GEN_KAPPA_KAZEIN ,
  genRj.GEN_BETA_GENOTIP,
  gp.ZGP_VREDNOST;


select * from zivali ziv where ziv.ZIV_ID_SEQ=2468263 ;



--vse rjave zivali na območju MCSELECT DISTINCT
select DISTINCT  ziv.ZIV_ID_SEQ, ziv.SIF_SPOL,  genRj.GEN_KAPPA_KAZEIN, genRj.GEN_BETA_GENOTIP
FROM
  GOVEDO.REJCI_MLEKARNA mc,
  govedo.zivali ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma,
  GOVEDO.GEN_RJAVA genRj
WHERE
genRj.GEN_ZIV_ID_SEQ=ziv.ZIV_ID_SEQ and
  ziv.CRE_SIFRA_CREDA       =mc.RM_SIF_LOKACIJA
AND mc.RM_REJEC_ODDAJA_MLEKO=1
and pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
and ziv.AKTIVNA=1
and pasma.PASMA_SPADA=1
;


select * from ALL_GEN_IND_27012017 gen, zivali ziv where gen.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ and ziv.SIF_SPOL=2;



--odbrane črede v tablei credegenotipizacija
truncate table credegenotipizacija;

drop table F4F_REJE_15032017;
SELECT DISTINCT
  govedo_splosne_procedure.PREBERI_NASLOV_REJCA_1VRSTA (crede.F4F_CRE_SIFRA_CREDA), crede.*
FROm
  F4F_REJE crede
;


--NAKNADNO IZBRANE REJE, 11.5.2017
--rejci, ki oddajajo mleko in so v AT kontroli in imajo rjave krave

CREATE Table reje_11052017 as
SELECT DISTINCT
  mc.RM_SIF_LOKACIJA creda,
  lok.*
FROM
  GOVEDO.REJCI_MLEKARNA mc,
  GOVEDO.LOKACIJE_STEVILO_ZIVALI lok,
  GOVEDO.STATUSI_ZIVALI_CREDA stat
  WHERE
  lok.LOKACIJA              =mc.RM_SIF_LOKACIJA
  and mc.RM_SIF_LOKACIJA=stat.CRE_SIFRA_CREDA
  and stat.SIF_STATUS_CREDA=5
  AND mc.RM_REJEC_ODDAJA_MLEKO=1
  and lok.DATUM_POROCILA='30.04.2017'
  and lok.RJAR not in 0
  and lok.LOKACIJA in (15575,2992,13087,3497,6149,5825,9624,7871,8976,4725,6670,8946,4257,7128,9621);
  
  --rejci, ki oddajajo mleko in so v AT kontroli (in imajo več kot 10 rajvih krav) --> živali - ali so v laktaciji ali ne
SELECT ziv.ZIV_ID_SEQ, ziv.IME_ZIVAL,ziv.SIF_BREJA ,ziv.CRE_SIFRA_CREDA,ziv.DRZ_ORIG_ZIVAL || ziv.STEV_ORIG_ZIVAL ID_zivali, 
ziv.DAT_ROJSTVO, ziv.ZIV_MATI_SEQ mati, ziv.ZIV_OCE_SEQ oce, ziv.DRZ_STEV_ORIG_MATI || ziv.ZIV_STEV_ORIG_MATI matiID, ziv.DRZ_STEV_ORIG_OCE || ziv.ZIV_STEV_ORIG_OCE oceID
,max(tel.ZAP_TELITEV) StTel, max(tel.DAT_TELITEV) ZadnjaTel, pv.VREDNOST_12_PV
FROM
Reje_11052017 reje, GOVEDO.ZIVALI ziv,
GOVEDO.ZIVALI_PASMA_SPADA pasma,
GOVEDO.TELITVE tel, GOVEDO.LAKTACIJE lak,
ARHIV.PLEMENSKE_VREDNOSTI pv
  WHERE
  ziv.ZIV_ID_SEQ=pv.PV_ZIV_ID_SEQ(+)
  and pv.SIFRA_LAST(+)=79
  and ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ 
  and tel.TEL_ID_SEQ=lak.LAK_TEL_ID_SEQ and
  ziv.CRE_SIFRA_CREDA=reje.CREDA
  and pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
  and ziv.SIF_SPOL=2
  and ziv.AKTIVNA=1
  and pasma.PASMA_SPADA=1
 -- and ziv.CRE_SIFRA_CREDA=32343
 --group by ziv.CRE_SIFRA_CREDA
--group by ziv.ZIV_ID_SEQ
group by ziv.ZIV_ID_SEQ, ziv.DRZ_ORIG_ZIVAL,ziv.STEV_ORIG_ZIVAL, ziv.DAT_ROJSTVO, ziv.ZIV_MATI_SEQ, ziv.ZIV_OCE_SEQ,
 ziv.DRZ_STEV_ORIG_MATI, ziv.ZIV_STEV_ORIG_MATI , ziv.DRZ_STEV_ORIG_OCE , ziv.ZIV_STEV_ORIG_OCE, ziv.sif_breja, ziv.CRE_SIFRA_CREDA, pv.VREDNOST_12_PV, ziv.IME_ZIVAL
;

--že genotipizirane živali
select ziv.DRZ_ORIG_ZIVAL || ziv.STEV_ORIG_ZIVAL from GENOTIPIZIRANE_ZIVALI gen, zivali ziv where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ;

--dodatno izbrane živali 11.5.2017
select distinct dod.CRE_SIFRA_CREDA, crede.KONTROLOR from DODATNA_ODBIRA_F4F dod, GOVEDO.CRE_OPISI crede where crede.CRE_SIFRA_LOK=dod.CRE_SIFRA_CREDA;


govedo_splosne_procedure.preberi_naslov_kont_1vrsta (ZGPT_KONTROLOR, NULL, 1, NULL, NULL) ime_kontrolor,
 govedo_splosne_procedure.preberi_naslov_kont_1vrsta ('3013', NULL, NULL, 1, NULL) ulica_kontrolorja,
 govedo_splosne_procedure.preberi_naslov_kont_1vrsta (ZGPT_KONTROLOR, NULL, NULL, NULL, 1) mesto_kontrolorja,