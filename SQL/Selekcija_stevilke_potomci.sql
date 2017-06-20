SELECT *
FROM GENOTIPIZIRANE_ZIVALI gen,
  GOVEDO.ZIVALI_PASMA_SPADA pasma
WHERE gen.ZIV_ID_SEQ  =pasma.ZIV_ID_SEQ
AND (pasma.PASMA_SPADA=2
OR pasma.pasma_spada  =222);
SELECT ziv.ZIV_ID_SEQ seq,
  ziv.DAT_ROJSTVO rojstvo ,
  ziv.SIF_SPOL spol,
  gp.ZGP_VREDNOST kappa,
  pasma.PASMA_SPADA pasma
FROM GOVEDO.ZIVALI_GP gp,
  GOVEDO.ZIVALI_PASMA_SPADA pasma,
  govedo.zivali ziv
WHERE ziv.ZIV_ID_SEQ =gp.ZGP_ID_SEQ
AND pasma.ZIV_ID_SEQ =gp.ZGP_ZIV_ID_SEQ
AND gp.ZGP_SGPL_SIFRA=100;
SELECT ziv.DRZ_ORIG_ZIVAL,
  tel.DAT_TELITEV,
  ose.OSE_ZIV_OCE_SEQ
FROM govedo.zivali ziv,
  GOVEDO.BIKOVSKE_MATERE bm,
  GOVEDO.TELITVE tel,
  GOVEDO.OSEMENITVE ose
WHERE ziv.ZIV_ID_SEQ                    = ose.OSE_ZIV_OCE_SEQ
AND ose.OSE_ID_SEQ                      =tel.TEL_OSE_ID_SEQ
AND bm.BM_ZIV_ID_SEQ                    =tel.TEL_ZIV_ID_SEQ
AND (extract(YEAR FROM tel.DAT_TELITEV) = 2016);
SELECT extract(YEAR FROM tel.DAT_TELITEV) FROM GOVEDO.TELITVE tel;
SELECT *
FROM govedo.zivali ziv,
  GOVEDO.OSEMENITVE ose,
  GOVEDO.BIKOVSKE_MATERE bm
WHERE ose.OSE_ZIV_ID_SEQ            =bm.BM_ID_SEQ
AND ziv.ZIV_ID_SEQ                  =ose.OSE_ZIV_OCE_SEQ
AND (extract(YEAR FROM ose.DAT_OSEM)=2016);
SELECT *
FROM GOVEDO.OSEMENITVE ose,
  GOVEDO.BIKOVSKE_MATERE bm
WHERE ose.OSE_ZIV_ID_SEQ            =bm.BM_ID_SEQ
AND (extract(YEAR FROM ose.DAT_OSEM)=2016);
SELECT ziv.DRZ_ORIG_ZIVAL,
  pasma.PASMA_SPADA,
  COUNT(*)
FROM GOVEDO.OSEMENITVE ose,
  GOVEDO.BIKOVSKE_MATERE bm,
  govedo.zivali ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma
WHERE ziv.ZIV_ID_SEQ                 =ose.OSE_ZIV_OCE_SEQ
AND bm.BM_ZIV_ID_SEQ                 =ose.OSE_ZIV_ID_SEQ
AND ( extract(YEAR FROM ose.DAT_OSEM)=2016 )
AND ziv.ZIV_ID_SEQ                   =pasma.ZIV_ID_SEQ
GROUP BY ziv.DRZ_ORIG_ZIVAL,
  pasma.PASMA_SPADA;
SELECT ziv.DRZ_ORIG_ZIVAL,
  par.ZA_PASMO_BM,
  COUNT(DISTINCT ziv.ZIV_ID_SEQ)
FROM govedo.zivali ziv,
  GOVEDO.PB_ZA_BM par
WHERE par.LETO    =2016
AND ziv.ZIV_ID_SEQ=par.PB_ZIV_ID_SEQ
GROUP BY ziv.DRZ_ORIG_ZIVAL,
  par.ZA_PASMO_BM;
--število moških in ženskih telet rojenih op letih, povprečje 2000-2016
SELECT stleta.pasma,
  stleta.spol,
  ROUND(AVG(stleta.stevilo),2) stevilo
FROM
  (SELECT extract(YEAR FROM tel.DAT_TELITEV),
    pasma.PASMA_SPADA pasma,
    ziv.SIF_SPOL spol,
    COUNT(tel.TEL_ID_SEQ) stevilo
  FROM GOVEDO.TELITVE tel,
    govedo.zivali ziv,
    GOVEDO.ZIVALI_PASMA_SPADA pasma
  WHERE pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
  AND ziv.ZIV_ID_SEQ    =tel.TEL_ZIV_ID_SEQ
  AND extract(YEAR FROM tel.DAT_TELITEV) BETWEEN 2000 AND 2016
  GROUP BY extract(YEAR FROM tel.DAT_TELITEV),
    ziv.SIF_SPOL,
    pasma.PASMA_SPADA
  ) stleta
GROUP BY stleta.pasma,
  stleta.spol;
--povprečje zadnjih 10let, po spolu
SELECT stleta.leto,
  stleta.pasma,
  stleta.spol,
  ROUND(AVG(stleta.stevilo),2) stevilo
FROM
  (SELECT extract(YEAR FROM ziv.DAT_ROJSTVO) leto,
    pasma.PASMA_SPADA pasma,
    ziv.SIF_SPOL spol,
    COUNT(ziv.ZIV_ID_SEQ) stevilo
  FROM govedo.zivali ziv,
    GOVEDO.ZIVALI_PASMA_SPADA pasma
  WHERE pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
  AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
  GROUP BY extract(YEAR FROM ziv.DAT_ROJSTVO),
    ziv.SIF_SPOL,
    pasma.PASMA_SPADA
  ) stleta
GROUP BY stleta.leto,
  stleta.pasma,
  stleta.spol;
--število po spolu in letih, vsako leto posebej
SELECT extract(YEAR FROM ziv.DAT_ROJSTVO),
  pasma.PASMA_SPADA pasma,
  ziv.SIF_SPOL spol,
  COUNT(ziv.ZIV_ID_SEQ) stevilo
FROM govedo.zivali ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma
WHERE pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
GROUP BY extract(YEAR FROM ziv.DAT_ROJSTVO),
  ziv.SIF_SPOL,
  pasma.PASMA_SPADA;
--koliko jih pride v laktacijo
SELECT extract(YEAR FROM ziv.DAT_ROJSTVO),
  pasma.PASMA_SPADA pasma,
  COUNT(ziv.ZIV_ID_SEQ) stevilo
FROM govedo.zivali ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma,
  GOVEDO.TELITVE tel,
  GOVEDO.LAKTACIJE lak
WHERE pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
AND ziv.SIF_SPOL  =2
AND ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ
AND tel.TEL_ID_SEQ=lak.LAK_TEL_ID_SEQ
GROUP BY extract(YEAR FROM ziv.DAT_ROJSTVO),
  pasma.PASMA_SPADA;
--združi v eno tabelo
SELECT st.leto,
  st.pasma,
  st.stevilo,
  lakt.stevilo stevilo_lakt,
  ROUND((lakt.stevilo/st.stevilo) ,2) procent_Lakt
FROM
  (SELECT extract(YEAR FROM ziv.DAT_ROJSTVO) leto,
    pasma.PASMA_SPADA pasma,
    COUNT(ziv.ZIV_ID_SEQ) stevilo
  FROM govedo.zivali ziv,
    GOVEDO.ZIVALI_PASMA_SPADA pasma
  WHERE pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
  AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
  AND ziv.SIF_SPOL=2
  GROUP BY extract(YEAR FROM ziv.DAT_ROJSTVO),
    pasma.PASMA_SPADA
  ) st
INNER JOIN
  (SELECT extract(YEAR FROM ziv.DAT_ROJSTVO) leto,
    pasma.PASMA_SPADA pasma,
    COUNT(DISTINCT ziv.ZIV_ID_SEQ) stevilo
  FROM govedo.zivali ziv,
    GOVEDO.ZIVALI_PASMA_SPADA pasma,
    GOVEDO.TELITVE tel,
    GOVEDO.LAKTACIJE lak
  WHERE pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
  AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
  AND ziv.SIF_SPOL  =2
  AND ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ
  AND tel.TEL_ID_SEQ=lak.LAK_TEL_ID_SEQ
  GROUP BY extract(YEAR FROM ziv.DAT_ROJSTVO),
    pasma.PASMA_SPADA
  ) lakt
ON st.pasma =lakt.pasma
AND st.leto =lakt.leto;
--telitve
SELECT extract(YEAR FROM ziv.DAT_ROJSTVO) leto,
  pasma.PASMA_SPADA pasma,
  ziv.SIF_SPOL spol,
  COUNT(DISTINCT ziv.ZIV_ID_SEQ) stevilo
FROM govedo.zivali ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma,
  GOVEDO.TELITVE tel
WHERE pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
AND ziv.SIF_SPOL  =2
AND ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ
GROUP BY extract(YEAR FROM ziv.DAT_ROJSTVO),
  pasma.PASMA_SPADA,
  ziv.sif_spol;
--bikovske matere po letih
SELECT
  --ziv.ziv_id_seq,
  extract(YEAR FROM ziv.DAT_ROJSTVO) letoRoj,
  -- extract(YEAR FROM bm.DAT_STATUS_BM) letoPriz,
  pasma.PASMA_SPADA pasma,
  COUNT(DISTINCT ziv.ziv_id_seq)
FROM govedo.zivali ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma,
  GOVEDO.BIKOVSKE_MATERE bm
WHERE pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
AND ziv.SIF_SPOL  =2
AND ziv.ZIV_ID_SEQ=bm.bm_ZIV_ID_SEQ
GROUP BY extract(YEAR FROM ziv.DAT_ROJSTVO),
  --extract(YEAR FROM bm.DAT_STATUS_BM),
  pasma.PASMA_SPADA,
  ziv.sif_spol;
--bikovske matere po  letih rojstva - kdaj so bile priznane
SELECT extract(YEAR FROM ziv.DAT_ROJSTVO) letoRoj,
  bikovM.letoPriz,
  COUNT(DISTINCT bikovM.SEQ)
FROM govedo.zivali ziv,
  (SELECT ZIV.ZIV_ID_SEQ SEQ,
    extract(YEAR FROM bm.DAT_STATUS_BM) letoPriz,
    pasma.PASMA_SPADA pasma
  FROM GOVEDO.BIKOVSKE_MATERE bm,
    zivali ziv,
    GOVEDO.ZIVALI_PASMA_SPADA pasma
  WHERE pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
  AND ziv.ZIV_ID_SEQ    =bm.BM_ZIV_ID_SEQ
  ) bikovM
WHERE bikovM.sEq=ziv.ziv_id_seq
GROUP BY extract(YEAR FROM ziv.DAT_ROJSTVO),
  bikovM.letoPriz ;
--plemenski biki po  letih
SELECT extract(YEAR FROM ziv.DAT_ROJSTVO) leto,
  pasma.PASMA_SPADA pasma,
  ziv.DRZ_ORIG_ZIVAL,
  pb.SIF_STAT_TEST_BIK,
  pb.SIF_UPORABA_PB,
  zs.ZSG_DOLGO_IME,
  COUNT(DISTINCT ziv.ZIV_ID_SEQ) stevilo_LetoRoj
FROM govedo.zivali ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma,
  GOVEDO.PLEMENSKI_BIKI pb,
  GOVEDO.SEZNAM_SIFRANTOV ss,
  GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs
WHERE zs.ZSG_SS_ID_SEQ=ss.SS_ID_SEQ
AND pb.SIF_UPORABA_PB =zs.ZSG_SIFRA
AND pasma.ZIV_ID_SEQ  =ziv.ZIV_ID_SEQ
AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
AND ziv.SIF_SPOL  =1
AND ziv.ZIV_ID_SEQ=pb.PB_ZIV_ID_SEQ
AND ss.SS_ID_SEQ  =17
GROUP BY extract(YEAR FROM ziv.DAT_ROJSTVO),
  pasma.PASMA_SPADA,
  pb.SIF_STAT_TEST_BIK,
  pb.SIF_UPORABA_PB,
  ziv.DRZ_ORIG_ZIVAL,
  zs.ZSG_DOLGO_IME;
SELECT ss.SS_KRATKO_IME,
  zs.ZSG_DOLGO_IME,
  zs.ZSG_SIFRA
FROM GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs,
  GOVEDO.SEZNAM_SIFRANTOV ss
WHERE ss.SS_ID_SEQ=zs.ZSG_SS_ID_SEQ
AND ss.SS_ID_SEQ  =17;
SELECT * FROM GOVEDO.PLEMENSKI_BIKI;
--mrtvorojeni
SELECT tl.TL_TELE_MRTVOROJENO,
  st.POMEN,
  tl.TL_SPOL,
  tl.TL_PASMA,
  extract (YEAR FROM tl.TL_DAT_ROJ),
  COUNT(DISTINCT tl.TL_ZIV_ID_SEQ)
FROM GOVEDO.TABELA_TETOVIRNI_LIST tl,
  GOVEDO.SIFRANT_STANJE_TELETA st
WHERE st.SIFRA_STANJE_TEL=tl.TL_TELE_MRTVOROJENO
AND extract(YEAR FROM tl.TL_DAT_ROJ) BETWEEN 2002 AND 2016
GROUP BY tl.TL_TELE_MRTVOROJENO,
  st.POMEN,
  tl.TL_SPOL,
  extract (YEAR FROM tl.TL_DAT_ROJ),
  tl.TL_PASMA;
--potomci BM in elitnih bikov
SELECT AVG( COUNT(DISTINCT ziv.ZIV_ID_SEQ))
FROM govedo.zivali ziv,
  GOVEDO.BIKOVSKE_MATERE bm,
  GOVEDO.PB_ZA_BM pb,
  GOVEDO.ZIVALI_PASMA_SPADA pasma
WHERE pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
AND ziv.ZIV_MATI_SEQ  = bm.BM_ZIV_ID_SEQ
AND ziv.ZIV_OCE_SEQ   =pb.PB_ZIV_ID_SEQ
AND pasma.PASMA_SPADA =1
AND ziv.SIF_SPOL      =1
AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2000 AND 2016
GROUP BY extract(YEAR FROM ziv.DAT_ROJSTVO),
  pasma.pasma_spada;
--potomci BM in elitnih bikov --> koliko po očetu
SELECT COUNT(DISTINCT ziv.ZIV_ID_SEQ),
  ziv.ZIV_OCE_SEQ,
  extract(YEAR FROM ziv.DAT_ROJSTVO)
FROM govedo.zivali ziv,
  GOVEDO.BIKOVSKE_MATERE bm,
  GOVEDO.PB_ZA_BM pb,
  GOVEDO.ZIVALI_PASMA_SPADA pasma
WHERE pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
AND ziv.ZIV_MATI_SEQ  = bm.BM_ZIV_ID_SEQ
AND ziv.ZIV_OCE_SEQ   =pb.PB_ZIV_ID_SEQ
AND pasma.PASMA_SPADA =1
AND ziv.SIF_SPOL      =1
AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2014 AND 2016
GROUP BY ziv.ZIV_OCE_SEQ,
  extract(YEAR FROM ziv.DAT_ROJSTVO);
--koliko pb za bm na leto
SELECT * FROM GOVEDO.PB_ZA_BM pb WHERE pb.LETO=2016;




--koliko potomcev po AI biku (oz. po vseh šifrah) - ne samo rjavi potomci
select round(avg(poLetu.stePoto)), poLetu.SIF_UPORABA_PB from
(select round(avg(poOce.stPot)) stePoto, poOce.SIF_UPORABA_PB, poOce.letoRoj from
(SELECT count(DISTINCT ziv.ZIV_ID_SEQ) stPot,
  extract(YEAR FROM ziv.DAT_ROJSTVO) letoRoj,
  pb.PB_ZIV_ID_SEQ,
  pb.SIF_UPORABA_PB
FROM govedo.zivali ziv,
  GOVEDO.PLEMENSKI_BIKI pb,
  GOVEDO.ZIVALI_PASMA_SPADA pasma
WHERE pb.PB_ZIV_ID_SEQ                  =ziv.ZIV_OCE_SEQ
AND extract(YEAR FROM ziv.DAT_ROJSTVO) between 2005 and 2015
and pb.PB_ZIV_ID_SEQ = pasma.ZIV_ID_SEQ
and pasma.PASMA_SPADA=1 group by   pb.PB_ZIV_ID_SEQ,
  pb.SIF_UPORABA_PB,  extract(YEAR FROM ziv.DAT_ROJSTVO)) poOce
  where poOce.stPot >=50  group by poOce.SIF_UPORABA_PB, poOce.letoRoj) poLetu group by poLetu.SIF_UPORABA_PB ;