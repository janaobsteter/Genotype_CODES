--trenutno stanje
SELECT stanje.KATEGORIJA_ZIVALI,
  COUNT(DISTINCT stanje.seq)
FROM GOVEDO.TOM_AKTIVNE_ZIVALI stanje,
  govedo.lokacije_stevilo_zivali lok
WHERE lok.LOKACIJA    =stanje.LOKACIJA
AND lok.VRSTA_KONTROLE='AP'
AND stanje.PASMA      =1
GROUP BY stanje.KATEGORIJA_ZIVALI;
--trenutno stanje - kakšne šifre imajo PB
SELECT pb.SIF_UPORABA_PB,
  COUNT(DISTINCT stanje.seq),
  ( extract(YEAR FROM stanje.TAZV_DATUM))
FROM GOVEDO.TOM_AKTIVNE_ZIVALI_VSE stanje,
  govedo.lokacije_stevilo_zivali lok,
  GOVEDO.PLEMENSKI_BIKI pb,
  GOVEDO.SEZNAM_SIFRANTOV ss,
  GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs
WHERE lok.LOKACIJA                         =stanje.LOKACIJA
AND lok.VRSTA_KONTROLE                     ='AP'
AND stanje.PASMA                           =1
AND stanje.SEQ                             =pb.PB_ZIV_ID_SEQ
AND stanje.KATEGORIJA_ZIVALI               ='PB'
AND zs.ZSG_SS_ID_SEQ                       =ss.SS_ID_SEQ
AND pb.SIF_UPORABA_PB                      =zs.ZSG_SIFRA
AND ss.SS_ID_SEQ                           =17
AND (extract(MONTH FROM stanje.TAZV_DATUM))='07'
AND ( extract(YEAR FROM stanje.TAZV_DATUM)) BETWEEN 2014 AND 2016
GROUP BY pb.SIF_UPORABA_PB,
  ( extract(YEAR FROM stanje.TAZV_DATUM)) ;
"""  --število ženskih živali v kontroliranih čredahSELECT  KATEGORIJA_ZIVALI,count (*) stFROM govedo.tom_aktivne_zivali_vse tz, govedo.lokacije_stevilo_zivali kWHERE tz.LOKACIJA=k.LOKACIJAAND k.datum_porocila= TO_DATE ('&3', 'ddmmyyyy')and TAZV_DATUM= '&datum_poročila'AND k.VRSTA_KONTROLE='AP'AND SPOL=2group by KATEGORIJA_ZIVALI;  """
SELECT stanje.KATEGORIJA_ZIVALI,
  stanje.seq
  --COUNT(DISTINCT stanje.seq),
  --extract(year from stanje.TAZV_DATUM)
FROM GOVEDO.TOM_AKTIVNE_ZIVALI_vse stanje,
  govedo.lokacije_stevilo_zivali lok
WHERE lok.LOKACIJA                         =stanje.LOKACIJA
AND lok.VRSTA_KONTROLE                     ='AP'
AND stanje.PASMA                           =1
AND (extract(MONTH FROM stanje.TAZV_DATUM))='07'
AND ( extract(YEAR FROM stanje.TAZV_DATUM)) BETWEEN 2014 AND 2016
GROUP BY stanje.KATEGORIJA_ZIVALI,
  extract(YEAR FROM stanje.TAZV_DATUM);
--aktivni biki po kategorijah, trenutno
SELECT COUNT(DISTINCT akt.SEQ),
  pb.SIF_UPORABA_PB
FROM GOVEDO.ZIVALI_PASMA_SPADA pasma ,
  GOVEDO.PLEMENSKI_BIKI pb,
  GOVEDO.SEZNAM_SIFRANTOV ss,
  GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs,
  GOVEDO.TOM_AKTIVNE_ZIVALI akt
WHERE akt.SEQ        =pb.PB_ZIV_ID_SEQ
AND akt.SEQ          =pasma.ZIV_ID_SEQ
AND pasma.PASMA_SPADA=1
AND zs.ZSG_SS_ID_SEQ =ss.SS_ID_SEQ
AND pb.SIF_UPORABA_PB=zs.ZSG_SIFRA
AND ss.SS_ID_SEQ     =17
GROUP BY pb.SIF_UPORABA_PB;
--aktivni biki po kategorijah, povprecje
SELECT COUNT(DISTINCT akt.SEQ),
  (extract(YEAR FROM akt.TAZV_DATUM)),
  pb.SIF_UPORABA_PB
FROM GOVEDO.ZIVALI_PASMA_SPADA pasma ,
  GOVEDO.PLEMENSKI_BIKI pb,
  GOVEDO.SEZNAM_SIFRANTOV ss,
  GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs,
  GOVEDO.TOM_AKTIVNE_ZIVALI_VSE akt
WHERE akt.SEQ                           =pb.PB_ZIV_ID_SEQ
AND akt.SEQ                             =pasma.ZIV_ID_SEQ
AND pasma.PASMA_SPADA                   =1
AND zs.ZSG_SS_ID_SEQ                    =ss.SS_ID_SEQ
AND pb.SIF_UPORABA_PB                   =zs.ZSG_SIFRA
AND ss.SS_ID_SEQ                        =17
AND (extract(MONTH FROM akt.TAZV_DATUM))='07'
AND ( extract(YEAR FROM akt.TAZV_DATUM)) BETWEEN 2014 AND 2016
GROUP BY pb.SIF_UPORABA_PB,
  (extract(YEAR FROM akt.TAZV_DATUM));
--koliko trenutno aktivnih ženskih živali potomk po doloceni sifri PB
SELECT COUNT(DISTINCT akt.SEQ),
  pb.SIF_UPORABA_PB
FROM GOVEDO.ZIVALI ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma ,
  GOVEDO.PLEMENSKI_BIKI pb,
  GOVEDO.SEZNAM_SIFRANTOV ss,
  GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs,
  GOVEDO.TOM_AKTIVNE_ZIVALI akt
WHERE akt.SEQ        =ziv.ZIV_ID_SEQ
AND ziv.ZIV_OCE_SEQ  =pb.PB_ZIV_ID_SEQ
AND akt.SEQ          =pasma.ZIV_ID_SEQ
AND pasma.PASMA_SPADA=1
AND zs.ZSG_SS_ID_SEQ =ss.SS_ID_SEQ
AND pb.SIF_UPORABA_PB=zs.ZSG_SIFRA
AND ss.SS_ID_SEQ     =17
GROUP BY pb.SIF_UPORABA_PB;
--koliko trenutno aktivnih ženskih živali potomk po doloceni sifri PB, povprečje skozi leta, (AP KONTROLA)
SELECT COUNT(DISTINCT akt.SEQ),
  pb.SIF_UPORABA_PB,
  (extract(YEAR FROM akt.TAZV_DATUM))
FROM GOVEDO.ZIVALI ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma ,
  GOVEDO.PLEMENSKI_BIKI pb,
  GOVEDO.SEZNAM_SIFRANTOV ss,
  GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs,
  GOVEDO.TOM_AKTIVNE_ZIVALI_VSE akt,
  GOVEDO.LOKACIJE_STEVILO_ZIVALI lok
WHERE akt.SEQ                           =ziv.ZIV_ID_SEQ
AND ziv.ZIV_OCE_SEQ                     =pb.PB_ZIV_ID_SEQ
AND akt.SEQ                             =pasma.ZIV_ID_SEQ
AND pasma.PASMA_SPADA                   =1
AND zs.ZSG_SS_ID_SEQ                    =ss.SS_ID_SEQ
AND (extract(MONTH FROM akt.TAZV_DATUM))='07'
AND ( extract(YEAR FROM akt.TAZV_DATUM)) BETWEEN 2014 AND 2016
AND pb.SIF_UPORABA_PB =zs.ZSG_SIFRA
AND ss.SS_ID_SEQ      =17
AND ziv.SIF_SPOL      =2
AND lok.LOKACIJA      =akt.LOKACIJA
AND lok.VRSTA_KONTROLE='AP'
GROUP BY pb.SIF_UPORABA_PB,
  (extract(YEAR FROM akt.TAZV_DATUM));
--koliko očetov trenutno aktivnih živali potomk po doloceni sifri PB, povprečje skozi leta, (AP KONTROLA)
SELECT COUNT(DISTINCT ziv.ZIV_OCE_SEQ),
  pb.SIF_UPORABA_PB,
  (extract(YEAR FROM akt.TAZV_DATUM))
FROM GOVEDO.ZIVALI ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma ,
  GOVEDO.PLEMENSKI_BIKI pb,
  GOVEDO.SEZNAM_SIFRANTOV ss,
  GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs,
  GOVEDO.TOM_AKTIVNE_ZIVALI_VSE akt,
  GOVEDO.LOKACIJE_STEVILO_ZIVALI lok
WHERE akt.SEQ                           =ziv.ZIV_ID_SEQ
AND ziv.ZIV_OCE_SEQ                     =pb.PB_ZIV_ID_SEQ
AND akt.SEQ                             =pasma.ZIV_ID_SEQ
AND pasma.PASMA_SPADA                   =1
AND zs.ZSG_SS_ID_SEQ                    =ss.SS_ID_SEQ
AND (extract(MONTH FROM akt.TAZV_DATUM))='07'
AND ( extract(YEAR FROM akt.TAZV_DATUM)) BETWEEN 2014 AND 2016
AND pb.SIF_UPORABA_PB=zs.ZSG_SIFRA
AND ss.SS_ID_SEQ     =17
  -- and ziv.SIF_SPOL=2
AND lok.LOKACIJA      =akt.LOKACIJA
AND lok.VRSTA_KONTROLE='AP'
GROUP BY pb.SIF_UPORABA_PB,
  (extract(YEAR FROM akt.TAZV_DATUM));
--koliko stare so živali
SELECT stanje.KATEGORIJA_ZIVALI,
  -- stanje.seq,
  COUNT(DISTINCT stanje.seq),
  ((extract(YEAR FROM stanje.TAZV_DATUM)) - (extract(YEAR FROM ziv.DAT_ROJSTVO))) starost
FROM GOVEDO.TOM_AKTIVNE_ZIVALI_vse stanje,
  govedo.lokacije_stevilo_zivali lok,
  govedo.zivali ziv
WHERE ziv.ZIV_ID_SEQ                       =stanje.seq
AND lok.LOKACIJA                           =stanje.LOKACIJA
AND lok.VRSTA_KONTROLE                     ='AP'
AND stanje.PASMA                           =1
AND (extract(MONTH FROM stanje.TAZV_DATUM))='07'
AND ( extract(YEAR FROM stanje.TAZV_DATUM)) BETWEEN 2010 AND 2016
GROUP BY ((extract(YEAR FROM stanje.TAZV_DATUM)) - (extract(YEAR FROM ziv.DAT_ROJSTVO))),
  stanje.KATEGORIJA_ZIVALI;
SELECT stanje.KATEGORIJA_ZIVALI,
  COUNT( DISTINCT stanje.seq),
  pb.SIF_STAT_TEST_BIK,
  pb.SIF_UPORABA_PB,
  (extract(MONTH FROM stanje.TAZV_DATUM)) mesec
  --COUNT(DISTINCT stanje.seq),
  --extract(year from stanje.TAZV_DATUM)
FROM GOVEDO.TOM_AKTIVNE_ZIVALI_vse stanje,
  govedo.lokacije_stevilo_zivali lok,
  GOVEDO.PLEMENSKI_BIKI pb
WHERE pb.PB_ZIV_ID_SEQ=stanje.seq
AND lok.LOKACIJA      =stanje.LOKACIJA
AND lok.VRSTA_KONTROLE='AP'
AND stanje.PASMA      =1
  --AND (extract(month from stanje.TAZV_DATUM))='07'
AND ( extract(YEAR FROM stanje.TAZV_DATUM)) = 2016
GROUP BY stanje.KATEGORIJA_ZIVALI,
  pb.SIF_STAT_TEST_BIK,
  pb.SIF_UPORABA_PB,
  (extract(MONTH FROM stanje.TAZV_DATUM)) ;
SELECT
SELECT * FROM GOVEDO.LOKACIJE_STEVILO_ZIVALI;
SELECT ziv.SIF_SPOL,
  COUNT(DISTINCT ziv.ZIV_ID_SEQ)
FROM govedo.zivali ziv
WHERE ( extract(YEAR FROM ziv.DAT_ROJSTVO) ) =2016
AND ziv.SP1_SIFRA_PASMA                      =1
GROUP BY ziv.SIF_SPOL;
SELECT ziv.SIF_SPOL,
  COUNT(DISTINCT ziv.ZIV_ID_SEQ)
FROM govedo.zivali ziv
WHERE ( extract(YEAR FROM ziv.DAT_ROJSTVO) ) =2015
AND ziv.SP1_SIFRA_PASMA                      =1
GROUP BY ziv.SIF_SPOL;
SELECT ziv.SIF_SPOL,
  COUNT(DISTINCT ziv.ZIV_ID_SEQ)
FROM govedo.zivali ziv
WHERE ( extract(YEAR FROM ziv.DAT_ROJSTVO) ) =2014
AND ziv.SP1_SIFRA_PASMA                      =1
GROUP BY ziv.SIF_SPOL;
SELECT ziv.SIF_SPOL,
  COUNT(DISTINCT ziv.ZIV_ID_SEQ)
FROM govedo.zivali ziv
WHERE ( extract(YEAR FROM ziv.DAT_ROJSTVO) ) =2010
AND ziv.SP1_SIFRA_PASMA                      =1
GROUP BY ziv.SIF_SPOL;
SELECT ziv.SIF_SPOL,
  COUNT(DISTINCT ziv.ZIV_ID_SEQ),
  (extract(YEAR FROM ziv.DAT_ROJSTVO))
FROM govedo.zivali ziv,
  GOVEDO.LOKACIJE_STEVILO_ZIVALI lok
WHERE ziv.CRE_SIFRA_CREDA=lok.LOKACIJA
AND (extract(YEAR FROM ziv.DAT_ROJSTVO)) BETWEEN 2007 AND 2016
AND ziv.SP1_SIFRA_PASMA=1
  --and lok.VRSTA_KONTROLE='AP'
GROUP BY ziv.SIF_SPOL,
  (extract(YEAR FROM ziv.DAT_ROJSTVO));
--trenutno stanje
SELECT stanje.KATEGORIJA_ZIVALI,
  akt_bm.sttel,
  --ziv.DAT_ROJSTVO,
  COUNT(akt_bm.BM_ZIV_ID_SEQ)
  --COUNT(DISTINCT stanje.seq)
FROM GOVEDO.TOM_AKTIVNE_ZIVALI stanje,
  (SELECT DISTINCT bm.BM_ZIV_ID_SEQ,
    MAX(t.ZAP_TELITEV) sttel
  FROM bikovske_matere bm,
    telitve t,
    zivali z
  WHERE bm.bm_ziv_id_seq               = z.ziv_id_seq
  AND bm.bm_ziv_id_seq                 = t.tel_ziv_id_seq
  AND bm.SIF_STATUS_BM                IN (1)
  AND z.sp1_sifra_pasma               IN ( 1)
  AND Z.AKTIVNA                        =1
  AND TO_CHAR (dat_status_bm, 'yyyy') IN ( 2016,2017)
  GROUP BY bm.BM_ZIV_ID_SEQ
  ) akt_bm
WHERE akt_bm.BM_ZIV_ID_SEQ =stanje.seq
  --AND bm.BM_ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
GROUP BY stanje.KATEGORIJA_ZIVALI,
  akt_bm.sttel ;
--trenutne bikovske matere
(
SELECT DISTINCT bm.BM_ZIV_ID_SEQ
FROM bikovske_matere bm,
  telitve t,
  zivali z
WHERE bm.bm_ziv_id_seq               = z.ziv_id_seq
AND bm.bm_ziv_id_seq                 = t.tel_ziv_id_seq
AND bm.SIF_STATUS_BM                IN (1)
AND z.sp1_sifra_pasma               IN ( 1)
AND Z.AKTIVNA                        =1
AND TO_CHAR (dat_status_bm, 'yyyy') IN ( 2016,2017)
) akt_bm;
SELECT DISTINCT *
FROM GOVEDO.PLEMENSKI_BIKI pb,
  GOVEDO.ZIVALI ziv
WHERE ziv.ZIV_ID_SEQ                      =pb.PB_ZIV_ID_SEQ
AND ((extract (YEAR FROM ziv.DAT_ROJSTVO))=2010)
AND pb.SIF_STAT_TEST_BIK                 IN (2);
--to so potomci načrtnih parjenj iz 2010, ki so bili pozitivno testirani
SELECT DISTINCT *
FROM GOVEDO.PLEMENSKI_BIKI pb,
  GOVEDO.ZIVALI ziv,
  (SELECT DISTINCT
    --COUNT(DISTINCT ziv.ZIV_ID_SEQ) st,
    ziv.ZIV_ID_SEQ,
    extract(YEAR FROM ziv.DAT_ROJSTVO) letoRoj,
    ziv.SIF_SPOL spol,
    pasma.PASMA_SPADA
  FROM govedo.zivali ziv,
    GOVEDO.OSEMENITVE ose,
    GOVEDO.BIKOVSKE_MATERE bm,
    GOVEDO.PB_ZA_BM pb,
    GOVEDO.ZIVALI_PASMA_SPADA pasma,
    GOVEDO.ODJAVA_PRIJAVA_ZIVALI sel
  WHERE ose.OSE_ZIV_ID_SEQ                 =ziv.ZIV_MATI_SEQ
  AND extract(YEAR FROM ose.DAT_OSEM)      = extract(YEAR FROM bm.DAT_STATUS_BM)
  AND bm.SIF_STATUS_BM                     =1
  AND pasma.ZIV_ID_SEQ                     =ziv.ZIV_ID_SEQ
  AND ziv.ZIV_MATI_SEQ                     = bm.BM_ZIV_ID_SEQ
  AND sel.OPZ_ZIV_ID_SEQ                   =ziv.ZIV_ID_SEQ
  AND (extract (YEAR FROM ziv.DAT_ROJSTVO))=2010
  AND ziv.ZIV_OCE_SEQ                      =pb.PB_ZIV_ID_SEQ
  AND ziv.SP1_SIFRA_PASMA                  =1
  AND ziv.SIF_SPOL                         =1
  ) potomci2010
  --GROUP BY  extract(YEAR FROM ziv.DAT_ROJSTVO),  ziv.sif_spol,  pasma.pasma_spada;
WHERE ziv.ZIV_ID_SEQ                      =pb.PB_ZIV_ID_SEQ
AND ziv.ZIV_ID_SEQ                        =potomci2010.ZIV_ID_SEQ
AND ((extract (YEAR FROM ziv.DAT_ROJSTVO))=2010)
AND pb.SIF_STAT_TEST_BIK                 IN (1);
;
--remont za krave
--trenutno stanje
SELECT stanje.KATEGORIJA_ZIVALI,
  COUNT(DISTINCT stanje.seq)
FROM GOVEDO.TOM_AKTIVNE_ZIVALI_vse stanje,
  (SELECT DISTINCT stanje.KATEGORIJA_ZIVALI,
    stanje.seq seq
  FROM GOVEDO.TOM_AKTIVNE_ZIVALI_vse stanje,
    govedo.lokacije_stevilo_zivali lok
  WHERE lok.LOKACIJA          =stanje.LOKACIJA
  AND lok.VRSTA_KONTROLE      ='AP'
  AND stanje.PASMA            =1
  AND stanje.KATEGORIJA_ZIVALI='K'
  AND stanje.TAZV_DATUM       ='31.08.2015'
  ) krave15,
  govedo.lokacije_stevilo_zivali lok
WHERE stanje.seq            = krave15.seq
AND lok.LOKACIJA            =stanje.LOKACIJA
AND lok.VRSTA_KONTROLE      ='AP'
AND stanje.PASMA            =1
AND stanje.KATEGORIJA_ZIVALI='K'
AND stanje.TAZV_DATUM       ='31.08.2016'
GROUP BY stanje.KATEGORIJA_ZIVALI;


--koliko hčera po enem plemenskem biku
select round(avg(poOcetih.hcere),1), pbi.SIF_UPORABA_PB from
(
SELECT ziv.ziv_oce_seq,
  pb.SIF_UPORABA_PB,
  COUNT (DISTINCT ziv.ziv_id_seq) hcere
FROM govedo.zivali ziv,
  GOVEDO.PLEMENSKI_BIKI pb
WHERE pb.PB_ZIV_ID_SEQ=ziv.ZIV_OCE_SEQ
GROUP BY ziv.ZIV_OCE_SEQ,
  pb.SIF_UPORABA_PB
) poOcetih, zivali ziva, GOVEDO.PLEMENSKI_BIKI pbi where ziva.ZIV_ID_SEQ=pbi.PB_ZIV_ID_SEQ and pbi.PB_ZIV_ID_SEQ=poOcetih.ziv_oce_seq and ziva.SP1_SIFRA_PASMA=1 group by pbi.SIF_UPORABA_PB;
