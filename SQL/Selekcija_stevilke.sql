select * from GENOTIPIZIRANE_ZIVALI gen, GOVEDO.ZIVALI_PASMA_SPADA pasma where gen.ZIV_ID_SEQ=pasma.ZIV_ID_SEQ and (pasma.PASMA_SPADA=2 or pasma.pasma_spada=222);

select ziv.ZIV_ID_SEQ seq, ziv.DAT_ROJSTVO rojstvo ,ziv.SIF_SPOL spol, gp.ZGP_VREDNOST kappa, pasma.PASMA_SPADA pasma from GOVEDO.ZIVALI_GP gp, GOVEDO.ZIVALI_PASMA_SPADA pasma, govedo.zivali ziv where ziv.ZIV_ID_SEQ=gp.ZGP_ID_SEQ and pasma.ZIV_ID_SEQ=gp.ZGP_ZIV_ID_SEQ and gp.ZGP_SGPL_SIFRA=100;

select ziv.DRZ_ORIG_ZIVAL, tel.DAT_TELITEV, ose.OSE_ZIV_OCE_SEQ from govedo.zivali ziv, GOVEDO.BIKOVSKE_MATERE bm, GOVEDO.TELITVE tel, GOVEDO.OSEMENITVE ose where ziv.ZIV_ID_SEQ= ose.OSE_ZIV_OCE_SEQ and ose.OSE_ID_SEQ=tel.TEL_OSE_ID_SEQ and bm.BM_ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ and (extract(year from tel.DAT_TELITEV) = 2016);


select extract(year from tel.DAT_TELITEV) from GOVEDO.TELITVE tel;


select * from govedo.zivali ziv, GOVEDO.OSEMENITVE ose, GOVEDO.BIKOVSKE_MATERE bm where ose.OSE_ZIV_ID_SEQ=bm.BM_ID_SEQ and ziv.ZIV_ID_SEQ=ose.OSE_ZIV_OCE_SEQ and (extract(year from ose.DAT_OSEM)=2016); 

select * from GOVEDO.OSEMENITVE ose, GOVEDO.BIKOVSKE_MATERE bm where ose.OSE_ZIV_ID_SEQ=bm.BM_ID_SEQ and (extract(year from ose.DAT_OSEM)=2016);

SELECT
  ziv.DRZ_ORIG_ZIVAL,
  pasma.PASMA_SPADA,
  COUNT(*)
FROM
  GOVEDO.OSEMENITVE ose,
  GOVEDO.BIKOVSKE_MATERE bm,
  govedo.zivali ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma
WHERE
  ziv.ZIV_ID_SEQ    =ose.OSE_ZIV_OCE_SEQ
AND bm.BM_ZIV_ID_SEQ=ose.OSE_ZIV_ID_SEQ
AND
  (
    extract(YEAR FROM ose.DAT_OSEM)=2016
  )
and ziv.ZIV_ID_SEQ=pasma.ZIV_ID_SEQ
GROUP BY
  ziv.DRZ_ORIG_ZIVAL,
  pasma.PASMA_SPADA;
  
  
  
select  ziv.DRZ_ORIG_ZIVAL, par.ZA_PASMO_BM, count(distinct ziv.ZIV_ID_SEQ) from govedo.zivali ziv, GOVEDO.PB_ZA_BM par where par.LETO=2016 and ziv.ZIV_ID_SEQ=par.PB_ZIV_ID_SEQ group by ziv.DRZ_ORIG_ZIVAL, par.ZA_PASMO_BM;

--število moških in ženskih telet rojenih op letih, povprečje 2000-2016
SELECT
  stleta.pasma,
  stleta.spol,
  round(AVG(stleta.stevilo),2) stevilo
FROM
  (
    SELECT
      extract(YEAR FROM tel.DAT_TELITEV),
      pasma.PASMA_SPADA pasma,
      ziv.SIF_SPOL spol,
      COUNT(DISTINCT tel.TEL_ID_SEQ) stevilo
    FROM
      GOVEDO.TELITVE tel,
      govedo.zivali ziv,
      GOVEDO.ZIVALI_PASMA_SPADA pasma,
      GOVEDO.LOKACIJE_STEVILO_ZIVALI lok
    WHERE
    lok.LOKACIJA=ziv.CRE_SIFRA_CREDA
    and lok.VRSTA_KONTROLE='AP'
    and  pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
    AND ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ
    AND extract(YEAR FROM tel.DAT_TELITEV) BETWEEN 2000 AND 2016
    GROUP BY
      extract(YEAR FROM tel.DAT_TELITEV),
      ziv.SIF_SPOL,
      pasma.PASMA_SPADA
  )
  stleta
GROUP BY
  stleta.pasma,
  stleta.spol;
  
  
--povprečje zadnjih 10let, po spolu
  SELECT
  stleta.pasma,
  stleta.spol,
  round(AVG(stleta.stevilo),2) stevilo_povp
FROM
  (
    SELECT
      extract(YEAR FROM ziv.DAT_ROJSTVO) leto,
      pasma.PASMA_SPADA pasma,
      ziv.SIF_SPOL spol,
      COUNT(distinct ziv.ZIV_ID_SEQ) stevilo
    FROM
      govedo.zivali ziv,
      GOVEDO.ZIVALI_PASMA_SPADA pasma,
      GOVEDO.LOKACIJE_STEVILO_ZIVALI LOK
    WHERE
    LOK.LOKACIJA=ziv.CRE_SIFRA_CREDA AND
    LOK.VRSTA_KONTROLE='AP' AND 
      pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
    AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2006 AND 2016
    GROUP BY
      extract(YEAR FROM ziv.DAT_ROJSTVO),
      ziv.SIF_SPOL,
      pasma.PASMA_SPADA
  )
  stleta
GROUP BY
  stleta.pasma,
  stleta.spol;

--število po spolu in letih, vsako leto posebej
    SELECT
      extract(YEAR FROM ziv.DAT_ROJSTVO),
      pasma.PASMA_SPADA pasma,
      ziv.SIF_SPOL spol,
      COUNT(ziv.ZIV_ID_SEQ) stevilo
    FROM
      govedo.zivali ziv,
      GOVEDO.ZIVALI_PASMA_SPADA pasma
    WHERE
      pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
    AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
    GROUP BY
      extract(YEAR FROM ziv.DAT_ROJSTVO),
      ziv.SIF_SPOL,
      pasma.PASMA_SPADA;
      
--koliko jih pride v laktacijo
 
  SELECT
      extract(YEAR FROM ziv.DAT_ROJSTVO),
      pasma.PASMA_SPADA pasma,
      COUNT(ziv.ZIV_ID_SEQ) stevilo
    FROM
      govedo.zivali ziv,
      GOVEDO.ZIVALI_PASMA_SPADA pasma,
      GOVEDO.TELITVE tel,
      GOVEDO.LAKTACIJE lak
    WHERE
      pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
    AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
    and ziv.SIF_SPOL=2
    and ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ
    and tel.TEL_ID_SEQ=lak.LAK_TEL_ID_SEQ
    GROUP BY
      extract(YEAR FROM ziv.DAT_ROJSTVO),
      pasma.PASMA_SPADA;
      

--združi v eno tabelo
select st.leto, st.pasma, st.stevilo, lakt.stevilo stevilo_lakt, round((lakt.stevilo/st.stevilo) ,2) procent_Lakt from
(    SELECT
      extract(YEAR FROM ziv.DAT_ROJSTVO) leto,
      pasma.PASMA_SPADA pasma,
      COUNT(ziv.ZIV_ID_SEQ) stevilo
    FROM
      govedo.zivali ziv,
      GOVEDO.ZIVALI_PASMA_SPADA pasma
    WHERE
      pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
    AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
    and ziv.SIF_SPOL=2
    GROUP BY
      extract(YEAR FROM ziv.DAT_ROJSTVO),
      pasma.PASMA_SPADA) st
inner join 
  (SELECT
      extract(YEAR FROM ziv.DAT_ROJSTVO) leto,
      pasma.PASMA_SPADA pasma,
      COUNT(distinct ziv.ZIV_ID_SEQ) stevilo
    FROM
      govedo.zivali ziv,
      GOVEDO.ZIVALI_PASMA_SPADA pasma,
      GOVEDO.TELITVE tel,
      GOVEDO.LAKTACIJE lak
    WHERE
      pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
    AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
    and ziv.SIF_SPOL=2
    and ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ
    and tel.TEL_ID_SEQ=lak.LAK_TEL_ID_SEQ
    GROUP BY
      extract(YEAR FROM ziv.DAT_ROJSTVO),
      pasma.PASMA_SPADA) lakt
    on st.pasma=lakt.pasma and st.leto =lakt.leto;
    
    
--telitve
SELECT
      extract(YEAR FROM ziv.DAT_ROJSTVO) leto,
      pasma.PASMA_SPADA pasma,
      ziv.SIF_SPOL spol,
      COUNT(distinct ziv.ZIV_ID_SEQ) stevilo
    FROM
      govedo.zivali ziv,
      GOVEDO.ZIVALI_PASMA_SPADA pasma,
      GOVEDO.TELITVE tel
    WHERE
      pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
    AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
    and ziv.SIF_SPOL=2
    and ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ
    GROUP BY
      extract(YEAR FROM ziv.DAT_ROJSTVO),
      pasma.PASMA_SPADA,
      ziv.sif_spol;
      
  --bikovske matere po letih
SELECT
  --ziv.ziv_id_seq,
  extract(YEAR FROM ziv.DAT_ROJSTVO) letoRoj,
 -- extract(YEAR FROM bm.DAT_STATUS_BM) letoPriz,
  pasma.PASMA_SPADA pasma,
  count(distinct ziv.ziv_id_seq)
FROM
  govedo.zivali ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma,
  GOVEDO.BIKOVSKE_MATERE bm
WHERE
  pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
AND ziv.SIF_SPOL  =2
AND ziv.ZIV_ID_SEQ=bm.bm_ZIV_ID_SEQ
GROUP BY extract(YEAR FROM ziv.DAT_ROJSTVO),
--extract(YEAR FROM bm.DAT_STATUS_BM),
pasma.PASMA_SPADA,
ziv.sif_spol;
  
  
  --bikovske matere po  letih rojstva - kdaj so bile priznane
select extract(year from ziv.DAT_ROJSTVO) letoRoj, bikovM.letoPriz, count(distinct bikovM.SEQ) 
from govedo.zivali ziv, 
(SELECT
ZIV.ZIV_ID_SEQ SEQ,
  extract(YEAR FROM bm.DAT_STATUS_BM) letoPriz,
  pasma.PASMA_SPADA pasma
FROM
  GOVEDO.BIKOVSKE_MATERE bm,
  zivali ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma
WHERE
  pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
AND ziv.ZIV_ID_SEQ=bm.BM_ZIV_ID_SEQ) bikovM
  where bikovM.sEq=ziv.ziv_id_seq
  group by extract(year from ziv.DAT_ROJSTVO), bikovM.letoPriz ;
      
  --plemenski biki po  letih
SELECT
      extract(YEAR FROM ziv.DAT_ROJSTVO) leto,
      pasma.PASMA_SPADA pasma,
      ziv.DRZ_ORIG_ZIVAL,
      pb.SIF_STAT_TEST_BIK,
      pb.SIF_UPORABA_PB,
      zs.ZSG_DOLGO_IME,
      COUNT(distinct ziv.ZIV_ID_SEQ) stevilo_LetoRoj
    FROM
      govedo.zivali ziv,
      GOVEDO.ZIVALI_PASMA_SPADA pasma,
      GOVEDO.PLEMENSKI_BIKI pb,
      GOVEDO.SEZNAM_SIFRANTOV ss,
      GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs
    WHERE
      zs.ZSG_SS_ID_SEQ=ss.SS_ID_SEQ and 
      pb.SIF_UPORABA_PB=zs.ZSG_SIFRA and
      pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
    AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
    and ziv.SIF_SPOL=1
    and ziv.ZIV_ID_SEQ=pb.PB_ZIV_ID_SEQ
    and ss.SS_ID_SEQ=17
    GROUP BY
      extract(YEAR FROM ziv.DAT_ROJSTVO),
      pasma.PASMA_SPADA,
      pb.SIF_STAT_TEST_BIK,
      pb.SIF_UPORABA_PB,
      ziv.DRZ_ORIG_ZIVAL,
      zs.ZSG_DOLGO_IME;
      
select ss.SS_KRATKO_IME, zs.ZSG_DOLGO_IME, zs.ZSG_SIFRA from GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs, GOVEDO.SEZNAM_SIFRANTOV ss where ss.SS_ID_SEQ=zs.ZSG_SS_ID_SEQ and ss.SS_ID_SEQ=17;
      

select * from GOVEDO.PLEMENSKI_BIKI;

--mrtvorojeni
SELECT
  tl.TL_TELE_MRTVOROJENO,
  st.POMEN,
  tl.TL_SPOL,
  tl.TL_PASMA,
  extract (year from tl.TL_DAT_ROJ),
  COUNT(DISTINCT tl.TL_ZIV_ID_SEQ)
FROM
  GOVEDO.TABELA_TETOVIRNI_LIST tl,
  GOVEDO.SIFRANT_STANJE_TELETA st
WHERE
  st.SIFRA_STANJE_TEL=tl.TL_TELE_MRTVOROJENO 
  and  extract(YEAR FROM tl.TL_DAT_ROJ) BETWEEN 2002 AND 2016
GROUP BY
   tl.TL_TELE_MRTVOROJENO,
  st.POMEN,
  tl.TL_SPOL,
  extract (year from tl.TL_DAT_ROJ),
  tl.TL_PASMA;


--potomci BM in elitnih bikov
SELECT
  COUNT(DISTINCT ziv.ZIV_ID_SEQ),
  extract(YEAR FROM ziv.DAT_ROJSTVO) letoRoj,
  ziv.SIF_SPOL spol,
  pasma.PASMA_SPADA
FROM
  govedo.zivali ziv,
  GOVEDO.BIKOVSKE_MATERE bm,
  GOVEDO.PB_ZA_BM pb,
  GOVEDO.ZIVALI_PASMA_SPADA pasma
WHERE
pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ and
  ziv.ZIV_MATI_SEQ = bm.BM_ZIV_ID_SEQ
AND ziv.ZIV_OCE_SEQ=pb.PB_ZIV_ID_SEQ
GROUP BY
  extract(YEAR FROM ziv.DAT_ROJSTVO),
  ziv.sif_spol,
  pasma.pasma_spada;
      
      
      
select count(distinct gen.ZIV_ID_SEQ) from JANAO.ALL_GEN_IND_27012017 gen, govedo.zivali ziv, GOVEDO.ZIVALI_PASMA_SPADA PASMA
where gen.ZIV_ID_SEQ=ziv.ziv_id_seq
AND ziv.ZIV_ID_SEQ=PASMA.ZIV_ID_SEQ
AND ziv.SIF_SPOL=1
AND PASMA.PASMA_SPADA=1
and (extract(year from ziv.DAT_ROJSTVO))=2015;