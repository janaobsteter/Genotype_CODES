                                                                                  --beta caseins
SELECT DISTINCT
  gen.GEN_ZIV_ID_SEQ ,
  gen.GEN_BCNAB,
  gen.GEN_BCN_A2,
  gen.GEN_BETA_GENOTIP,
  ziv.SIF_SPOL
FROM
  GOVEDO.GEN_RJAVA gen,
  zivali ziv
WHERE
  ziv.ZIV_ID_SEQ=gen.GEN_ZIV_ID_SEQ;
  

select gen.ZIV_ID_SEQ, ziv.DRZ_ORIG_ZIVAL || ziv.STEV_ORIG_ZIVAL from GENOTIPIZIRANE_ZIVALI gen, zivali ziv where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ;
select * from GENOTIPIZIRANE_ZIVALI gen where gen.ZIV_ID_SEQ=3786992;


SELECT DISTINCT
  ziv.DRZ_ORIG_ZIVAL
  || ziv.STEV_ORIG_ZIVAL,
  genziv.ZIV_ID_SEQ,
  genziv.GEN_DATUM,
  genziv.GEN_CHIP
FROM
  GOVEDO.GEN_RJAVA gen,
  GENOTIPIZIRANE_ZIVALI genziv,
  govedo.zivali ziv
WHERE
  ziv.ZIV_ID_SEQ      =gen.GEN_ZIV_ID_SEQ
AND gen.GEN_BCNAB    IS NULL
AND gen.GEN_ZIV_ID_SEQ=genziv.ZIV_ID_SEQ;
SELECT
  *
FROM
  GOVEDO.zivali_gp gen;
SELECT
  *
FROM
  GOVEDO.SIFRANT_GP_LASTNOSTI;
  
  --kapa = 100, beta = 103
SELECT
  ziv.ZIV_ID_SEQ seq,
  ziv.DRZ_ORIG_ZIVAL || ziv.ZIV_ID_SEQ,
  extract(year from ziv.DAT_ROJSTVO) rojstvo,
  ziv.SIF_SPOL spol,
  gp.ZGP_VREDNOST kappa,
  pasma.PASMA_SPADA pasma
FROM
  govedo.zivali_gp gp,
  govedo.zivali ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma
WHERE
  pasma.ZIV_ID_SEQ    =ziv.ZIV_ID_SEQ
AND gp.ZGP_SGPL_SIFRA =100
AND ziv.ZIV_ID_SEQ    =gp.ZGP_ZIV_ID_SEQ
and pasma.pasma_spada=1;