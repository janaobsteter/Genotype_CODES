
--trenutno stanje
SELECT
  stanje.KATEGORIJA_ZIVALI,
  COUNT(DISTINCT stanje.seq)
FROM
  GOVEDO.TOM_AKTIVNE_ZIVALI stanje,
   govedo.lokacije_stevilo_zivali lok
WHERE
lok.LOKACIJA=stanje.LOKACIJA
and lok.VRSTA_KONTROLE='AP'
 and stanje.PASMA=1
GROUP BY
  stanje.KATEGORIJA_ZIVALI;
  
  
  --trenutno stanje - kakšne šifre imajo PB
SELECT
  pb.SIF_UPORABA_PB,
  COUNT(DISTINCT stanje.seq),
   (  extract(year from stanje.TAZV_DATUM)) 
FROM
  GOVEDO.TOM_AKTIVNE_ZIVALI_VSE stanje,
   govedo.lokacije_stevilo_zivali lok, GOVEDO.PLEMENSKI_BIKI pb,
      GOVEDO.SEZNAM_SIFRANTOV ss,
      GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs
WHERE
lok.LOKACIJA=stanje.LOKACIJA
and lok.VRSTA_KONTROLE='AP'
 and stanje.PASMA=1
 and stanje.SEQ=pb.PB_ZIV_ID_SEQ
 and stanje.KATEGORIJA_ZIVALI='PB'
 and       zs.ZSG_SS_ID_SEQ=ss.SS_ID_SEQ and 
      pb.SIF_UPORABA_PB=zs.ZSG_SIFRA
    AND  ss.SS_ID_SEQ=17
     AND (extract(month from stanje.TAZV_DATUM))='07'
 and (  extract(year from stanje.TAZV_DATUM)) between 2014 and 2016
GROUP BY
pb.SIF_UPORABA_PB,  (  extract(year from stanje.TAZV_DATUM)) ;
  
"""  
--število ženskih živali v kontroliranih čredah
SELECT  KATEGORIJA_ZIVALI,
count (*) st
FROM govedo.tom_aktivne_zivali_vse tz, govedo.lokacije_stevilo_zivali k
WHERE 
tz.LOKACIJA=k.LOKACIJA
AND k.datum_porocila= TO_DATE ('&3', 'ddmmyyyy')
and TAZV_DATUM= '&datum_poročila'
AND k.VRSTA_KONTROLE='AP'
AND SPOL=2
group by KATEGORIJA_ZIVALI;  
""" 

SELECT
  stanje.KATEGORIJA_ZIVALI,
  stanje.seq
  --COUNT(DISTINCT stanje.seq),
  --extract(year from stanje.TAZV_DATUM)
FROM
  GOVEDO.TOM_AKTIVNE_ZIVALI_vse stanje,
   govedo.lokacije_stevilo_zivali lok
WHERE
lok.LOKACIJA=stanje.LOKACIJA
and lok.VRSTA_KONTROLE='AP'
 and stanje.PASMA=1
 AND (extract(month from stanje.TAZV_DATUM))='07'
 and (  extract(year from stanje.TAZV_DATUM)) between 2014 and 2016
GROUP BY
  stanje.KATEGORIJA_ZIVALI,
  extract(year from stanje.TAZV_DATUM);


--aktivni biki po kategorijah, trenutno
select  count(distinct akt.SEQ),  pb.SIF_UPORABA_PB from  GOVEDO.ZIVALI_PASMA_SPADA pasma , GOVEDO.PLEMENSKI_BIKI pb,
      GOVEDO.SEZNAM_SIFRANTOV ss,
      GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs, GOVEDO.TOM_AKTIVNE_ZIVALI akt
where  akt.SEQ=pb.PB_ZIV_ID_SEQ and akt.SEQ=pasma.ZIV_ID_SEQ and pasma.PASMA_SPADA=1 
and       zs.ZSG_SS_ID_SEQ=ss.SS_ID_SEQ and 
      pb.SIF_UPORABA_PB=zs.ZSG_SIFRA
    AND  ss.SS_ID_SEQ=17

    group by pb.SIF_UPORABA_PB;
    

--aktivni biki po kategorijah, povprecje
select  count(distinct akt.SEQ),  (extract(year from akt.TAZV_DATUM)), pb.SIF_UPORABA_PB from  GOVEDO.ZIVALI_PASMA_SPADA pasma , GOVEDO.PLEMENSKI_BIKI pb,
      GOVEDO.SEZNAM_SIFRANTOV ss, 
      GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs,
      GOVEDO.TOM_AKTIVNE_ZIVALI_VSE akt
where  akt.SEQ=pb.PB_ZIV_ID_SEQ and akt.SEQ=pasma.ZIV_ID_SEQ and pasma.PASMA_SPADA=1 
and       zs.ZSG_SS_ID_SEQ=ss.SS_ID_SEQ and 
      pb.SIF_UPORABA_PB=zs.ZSG_SIFRA
    AND  ss.SS_ID_SEQ=17
     AND (extract(month from akt.TAZV_DATUM))='07'
 and (  extract(year from akt.TAZV_DATUM)) between 2014 and 2016
    group by pb.SIF_UPORABA_PB, (extract(year from akt.TAZV_DATUM));
    
    
--koliko trenutno aktivnih ženskih živali potomk po doloceni sifri PB
select  count(distinct akt.SEQ),  pb.SIF_UPORABA_PB from  GOVEDO.ZIVALI ziv, GOVEDO.ZIVALI_PASMA_SPADA pasma , GOVEDO.PLEMENSKI_BIKI pb,
      GOVEDO.SEZNAM_SIFRANTOV ss,
      GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs, GOVEDO.TOM_AKTIVNE_ZIVALI akt
where akt.SEQ=ziv.ZIV_ID_SEQ and ziv.ZIV_OCE_SEQ=pb.PB_ZIV_ID_SEQ and akt.SEQ=pasma.ZIV_ID_SEQ and pasma.PASMA_SPADA=1 
and       zs.ZSG_SS_ID_SEQ=ss.SS_ID_SEQ and 
      pb.SIF_UPORABA_PB=zs.ZSG_SIFRA
    AND  ss.SS_ID_SEQ=17
    group by pb.SIF_UPORABA_PB;
    
        
--koliko trenutno aktivnih ženskih živali potomk po doloceni sifri PB, povprečje skozi leta, (AP KONTROLA)
select  count(distinct akt.SEQ),  pb.SIF_UPORABA_PB, (extract(year from akt.TAZV_DATUM)) from  GOVEDO.ZIVALI ziv, GOVEDO.ZIVALI_PASMA_SPADA pasma , GOVEDO.PLEMENSKI_BIKI pb,
      GOVEDO.SEZNAM_SIFRANTOV ss, 
      GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs, GOVEDO.TOM_AKTIVNE_ZIVALI_VSE akt, GOVEDO.LOKACIJE_STEVILO_ZIVALI lok
where akt.SEQ=ziv.ZIV_ID_SEQ and ziv.ZIV_OCE_SEQ=pb.PB_ZIV_ID_SEQ and akt.SEQ=pasma.ZIV_ID_SEQ and pasma.PASMA_SPADA=1 
and       zs.ZSG_SS_ID_SEQ=ss.SS_ID_SEQ 
     AND (extract(month from akt.TAZV_DATUM))='07'
 and (  extract(year from akt.TAZV_DATUM)) between 2014 and 2016
    and  pb.SIF_UPORABA_PB=zs.ZSG_SIFRA
    AND  ss.SS_ID_SEQ=17
    and ziv.SIF_SPOL=2
    and lok.LOKACIJA=akt.LOKACIJA
    and lok.VRSTA_KONTROLE='AP'
    group by pb.SIF_UPORABA_PB, (extract(year from akt.TAZV_DATUM));


--koliko očetov trenutno aktivnih živali potomk po doloceni sifri PB, povprečje skozi leta, (AP KONTROLA)
select  count(distinct ziv.ZIV_OCE_SEQ),  pb.SIF_UPORABA_PB, (extract(year from akt.TAZV_DATUM)) from  GOVEDO.ZIVALI ziv, GOVEDO.ZIVALI_PASMA_SPADA pasma , GOVEDO.PLEMENSKI_BIKI pb,
      GOVEDO.SEZNAM_SIFRANTOV ss, 
      GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs, GOVEDO.TOM_AKTIVNE_ZIVALI_VSE akt, GOVEDO.LOKACIJE_STEVILO_ZIVALI lok
where akt.SEQ=ziv.ZIV_ID_SEQ and ziv.ZIV_OCE_SEQ=pb.PB_ZIV_ID_SEQ and akt.SEQ=pasma.ZIV_ID_SEQ and pasma.PASMA_SPADA=1 
and       zs.ZSG_SS_ID_SEQ=ss.SS_ID_SEQ 
     AND (extract(month from akt.TAZV_DATUM))='07'
 and (  extract(year from akt.TAZV_DATUM)) between 2014 and 2016
    and  pb.SIF_UPORABA_PB=zs.ZSG_SIFRA
    AND  ss.SS_ID_SEQ=17
   -- and ziv.SIF_SPOL=2
    and lok.LOKACIJA=akt.LOKACIJA
    and lok.VRSTA_KONTROLE='AP'
    group by pb.SIF_UPORABA_PB, (extract(year from akt.TAZV_DATUM));

--koliko stare so živali
SELECT
  stanje.KATEGORIJA_ZIVALI,
 -- stanje.seq,
  COUNT(DISTINCT stanje.seq),
((extract(year from stanje.TAZV_DATUM)) - (extract(year from ziv.DAT_ROJSTVO))) starost
FROM
  GOVEDO.TOM_AKTIVNE_ZIVALI_vse stanje,
   govedo.lokacije_stevilo_zivali lok,
   govedo.zivali ziv
WHERE
ziv.ZIV_ID_SEQ=stanje.seq and
lok.LOKACIJA=stanje.LOKACIJA
and lok.VRSTA_KONTROLE='AP'
 and stanje.PASMA=1
 AND (extract(month from stanje.TAZV_DATUM))='07'
 and (  extract(year from stanje.TAZV_DATUM)) between 2010 and 2016
GROUP BY
((extract(year from stanje.TAZV_DATUM)) - (extract(year from ziv.DAT_ROJSTVO))),
  stanje.KATEGORIJA_ZIVALI;
  

SELECT
  stanje.KATEGORIJA_ZIVALI,
  count( distinct stanje.seq),
  pb.SIF_STAT_TEST_BIK,
  pb.SIF_UPORABA_PB,
  (extract(month from stanje.TAZV_DATUM)) mesec
  --COUNT(DISTINCT stanje.seq),
  --extract(year from stanje.TAZV_DATUM)
FROM
  GOVEDO.TOM_AKTIVNE_ZIVALI_vse stanje,
   govedo.lokacije_stevilo_zivali lok,
   GOVEDO.PLEMENSKI_BIKI pb
WHERE
pb.PB_ZIV_ID_SEQ=stanje.seq
and lok.LOKACIJA=stanje.LOKACIJA
and lok.VRSTA_KONTROLE='AP'
 and stanje.PASMA=1
 --AND (extract(month from stanje.TAZV_DATUM))='07'
 and (  extract(year from stanje.TAZV_DATUM)) = 2016
GROUP BY
  stanje.KATEGORIJA_ZIVALI,
  pb.SIF_STAT_TEST_BIK,
    pb.SIF_UPORABA_PB,
    (extract(month from stanje.TAZV_DATUM))
;

select 
  
  
  
  
  
SELECT
  *
FROM
  GOVEDO.LOKACIJE_STEVILO_ZIVALI;
SELECT
  ziv.SIF_SPOL,
  COUNT(DISTINCT ziv.ZIV_ID_SEQ)
FROM
  govedo.zivali ziv
WHERE
  (
    extract(YEAR FROM ziv.DAT_ROJSTVO)
  )
                       =2016
AND ziv.SP1_SIFRA_PASMA=1
GROUP BY
  ziv.SIF_SPOL;
SELECT
  ziv.SIF_SPOL,
  COUNT(DISTINCT ziv.ZIV_ID_SEQ)
FROM
  govedo.zivali ziv
WHERE
  (
    extract(YEAR FROM ziv.DAT_ROJSTVO)
  )
                       =2015
AND ziv.SP1_SIFRA_PASMA=1
GROUP BY
  ziv.SIF_SPOL;
SELECT
  ziv.SIF_SPOL,
  COUNT(DISTINCT ziv.ZIV_ID_SEQ)
FROM
  govedo.zivali ziv
WHERE
  (
    extract(YEAR FROM ziv.DAT_ROJSTVO)
  )
                       =2014
AND ziv.SP1_SIFRA_PASMA=1
GROUP BY
  ziv.SIF_SPOL;
SELECT
  ziv.SIF_SPOL,
  COUNT(DISTINCT ziv.ZIV_ID_SEQ)
FROM
  govedo.zivali ziv
WHERE
  (
    extract(YEAR FROM ziv.DAT_ROJSTVO)
  )
                       =2010
AND ziv.SP1_SIFRA_PASMA=1
GROUP BY
  ziv.SIF_SPOL;
  
SELECT  ziv.SIF_SPOL,
  COUNT(DISTINCT ziv.ZIV_ID_SEQ),(extract(YEAR FROM ziv.DAT_ROJSTVO))
FROM
  govedo.zivali ziv,
  GOVEDO.LOKACIJE_STEVILO_ZIVALI lok
WHERE
ziv.CRE_SIFRA_CREDA=lok.LOKACIJA
and  (extract(YEAR FROM ziv.DAT_ROJSTVO)) between 2007 and 2016
AND ziv.SP1_SIFRA_PASMA=1
--and lok.VRSTA_KONTROLE='AP'
GROUP BY
  ziv.SIF_SPOL,
  (extract(YEAR FROM ziv.DAT_ROJSTVO));
  
  



--trenutno stanje
SELECT
  stanje.KATEGORIJA_ZIVALI,
  akt_bm.sttel,
  --ziv.DAT_ROJSTVO,
 count(akt_bm.BM_ZIV_ID_SEQ)
  --COUNT(DISTINCT stanje.seq)
FROM
  GOVEDO.TOM_AKTIVNE_ZIVALI stanje,
(SELECT
  distinct bm.BM_ZIV_ID_SEQ, max(t.ZAP_TELITEV) sttel
FROM
  bikovske_matere bm,
  telitve t,
  zivali z
WHERE
  bm.bm_ziv_id_seq                   = z.ziv_id_seq
AND bm.bm_ziv_id_seq                 = t.tel_ziv_id_seq
AND bm.SIF_STATUS_BM                IN (1)
AND z.sp1_sifra_pasma               IN ( 1)
AND Z.AKTIVNA                        =1
AND TO_CHAR (dat_status_bm, 'yyyy') IN ( 2016,2017) group by bm.BM_ZIV_ID_SEQ) akt_bm
WHERE
  akt_bm.BM_ZIV_ID_SEQ  =stanje.seq
--AND bm.BM_ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
  GROUP BY stanje.KATEGORIJA_ZIVALI, akt_bm.sttel
  ;


--trenutne bikovske matere
(SELECT
  distinct bm.BM_ZIV_ID_SEQ
FROM
  bikovske_matere bm,
  telitve t,
  zivali z
WHERE
  bm.bm_ziv_id_seq                   = z.ziv_id_seq
AND bm.bm_ziv_id_seq                 = t.tel_ziv_id_seq
AND bm.SIF_STATUS_BM                IN (1)
AND z.sp1_sifra_pasma               IN ( 1)
AND Z.AKTIVNA                        =1
AND TO_CHAR (dat_status_bm, 'yyyy') IN ( 2016,2017)) akt_bm;


select distinct * from GOVEDO.PLEMENSKI_BIKI pb, GOVEDO.ZIVALI ziv where ziv.ZIV_ID_SEQ=pb.PB_ZIV_ID_SEQ and ((extract (year from ziv.DAT_ROJSTVO))=2010) and pb.SIF_STAT_TEST_BIK in (2);


--to so potomci načrtnih parjenj iz 2010, ki so bili pozitivno testirani
select distinct * from GOVEDO.PLEMENSKI_BIKI pb, GOVEDO.ZIVALI ziv, 
(SELECT distinct
  --COUNT(DISTINCT ziv.ZIV_ID_SEQ) st,
  ziv.ZIV_ID_SEQ,
  extract(YEAR FROM ziv.DAT_ROJSTVO) letoRoj,
  ziv.SIF_SPOL spol,
  pasma.PASMA_SPADA
FROM
  govedo.zivali ziv,GOVEDO.OSEMENITVE ose,
  GOVEDO.BIKOVSKE_MATERE bm,
  GOVEDO.PB_ZA_BM pb,
  GOVEDO.ZIVALI_PASMA_SPADA pasma,
  GOVEDO.ODJAVA_PRIJAVA_ZIVALI sel
WHERe
ose.OSE_ZIV_ID_SEQ=ziv.ZIV_MATI_SEQ and 
extract(year from ose.DAT_OSEM) = extract(year from bm.DAT_STATUS_BM)
and bm.SIF_STATUS_BM=1 and
pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ and
  ziv.ZIV_MATI_SEQ = bm.BM_ZIV_ID_SEQ
  and sel.OPZ_ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
  and (extract (year from ziv.DAT_ROJSTVO))=2010
AND ziv.ZIV_OCE_SEQ=pb.PB_ZIV_ID_SEQ
and ziv.SP1_SIFRA_PASMA=1
and ziv.SIF_SPOL=1) potomci2010
--GROUP BY  extract(YEAR FROM ziv.DAT_ROJSTVO),  ziv.sif_spol,  pasma.pasma_spada;
where ziv.ZIV_ID_SEQ=pb.PB_ZIV_ID_SEQ 
and ziv.ZIV_ID_SEQ=potomci2010.ZIV_ID_SEQ
and ((extract (year from ziv.DAT_ROJSTVO))=2010) 
and pb.SIF_STAT_TEST_BIK in (1);
;


--remont za krave

--trenutno stanje

SELECT
  stanje.KATEGORIJA_ZIVALI,
  count(DISTINCT stanje.seq)
FROM
  GOVEDO.TOM_AKTIVNE_ZIVALI_vse stanje,
   
(SELECT
  distinct stanje.KATEGORIJA_ZIVALI,
   stanje.seq seq
FROM
  GOVEDO.TOM_AKTIVNE_ZIVALI_vse stanje,
   govedo.lokacije_stevilo_zivali lok
WHERE
lok.LOKACIJA=stanje.LOKACIJA
and lok.VRSTA_KONTROLE='AP'
 and stanje.PASMA=1
 and stanje.KATEGORIJA_ZIVALI='K'
 and stanje.TAZV_DATUM='31.08.2015'
) krave15,
  govedo.lokacije_stevilo_zivali lok
WHERE
stanje.seq = krave15.seq
and lok.LOKACIJA=stanje.LOKACIJA
and lok.VRSTA_KONTROLE='AP'
 and stanje.PASMA=1
 and stanje.KATEGORIJA_ZIVALI='K'
 and stanje.TAZV_DATUM='31.08.2016'
GROUP BY
  stanje.KATEGORIJA_ZIVALI;
  
  
  

