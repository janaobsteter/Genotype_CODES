DROP TABLE OCETJE;
CREATE
  TABLE ocetje AS
SELECT DISTINCT
  gen.ZIV_ID_SEQ,
  ziv.ZIV_OCE_SEQ,
  ziv.ZIV_MATI_SEQ
FROM
  GENOTIPIZIRANE_ZIVALI gen,
  GOVEDO.ZIVALI ziv
WHERE
  gen.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ; --ocetje genotipiziranih zivali,ki so tudi
-- genotipizirani
CREATE
  TABLE GENOCETJE AS
SELECT DISTINCT
  ziv.ZIV_ID_SEQ,
  ziv.ZIV_MATI_SEQ,
  ziv.ZIV_OCE_SEQ
FROM
  OCETJE,
  PARENTAL_SNP800 gen,
  govedo.zivali ziv
WHERE
  OCETJE.ZIV_OCE_SEQ =gen.ZGP_ZIV_ID_SEQ
AND OCETJE.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ;

--pridobi očete genotipiziranih živali, ki imajo genotipizirano mamo
SELECT DISTINCT
  ziv.ZIV_ID_SEQ,
  ziv.ZIV_OCE_SEQ,
  OCETJE.ZIV_MATI_SEQ
FROM
  OCETJE,
  govedo.zivali ziv
WHERE
  OCETJE.ZIV_OCE_SEQ=ziv.ZIV_ID_SEQ;
SELECT DISTINCT
  ocetje.ZIV_ID_SEQ,
  ocetje.ZIV_OCE_SEQ,
  ziv.ZIV_ID_SEQ,
  ziv.DRZ_ORIG_ZIVAL
  || ziv.STEV_ORIG_ZIVAL
FROM
  (
    SELECT
      ziv.ZIV_ID_SEQ,
      ziv.ZIV_OCE_SEQ,
      ziv.ZIV_MATI_SEQ,
      ziv.DRZ_STEV_ORIG_OCE
    FROM
      govedo.zivali ziv
    WHERE
      ziv.ZIV_ID_SEQ=ANY(3622454, 4347149, 4460003, 4155287, 4434527, 4211966,
      4077913, 3881061, 3963929, 4148526, 4212593, 4306342, 4434528, 4232175,
      4278029, 4296189, 4355711)
  )
  ocetje,
  GOVEDO.ZIVALI ziv
WHERE
  ziv.ZIV_ID_SEQ=ocetje.ZIV_OCE_SEQ;
  
  
  
SELECT
  *
FROM
  govedo.zivali ziv
WHERE
  ziv.ZIV_ID_SEQ=ANY(4212593, 4306342, 3740925);
  
  
  
  
SELECT DISTINCT
  ocet.ZIV_ID_SEQ,
  ocet.ID_zivali,
  ocet.ZIV_OCE_SEQ,
  ziv.DRZ_ORIG_ZIVAL
  || ziv.STEV_ORIG_ZIVAL AS id_oceta,
  geno.GEN_CHIP
FROM
  (
    SELECT DISTINCT
      gen.ZIV_ID_SEQ,
      ziv.DRZ_ORIG_ZIVAL
      || ziv.STEV_ORIG_ZIVAL AS ID_zivali,
      ziv.ZIV_OCE_SEQ,
      gen.GEN_CHIP
    FROM
      govedo.zivali ziv,
      GENOTIPIZIRANE_ZIVALI gen
    WHERE
      gen.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
  )
  ocet,
  govedo.zivali ziv,
  GENOTIPIZIRANE_ZIVALI geno
WHERE
  ocet.ZIV_OCE_SEQ  =ziv.ZIV_ID_SEQ
AND ocet.ZIV_OCE_SEQ=geno.ZIV_ID_SEQ(+)
AND geno.GEN_CHIP  IS NULL;


select * from GOVEDO.ZIVALI ziv where ziv.STEV_ORIG_ZIVAL=1029117;

, 14234240, 23903151, 362960,
  43745627, 43867222, 44038436, 53575045, 54333510, 559287, 559296, 63616529,
  74234237
  
  
--find grandfather
 select ziv.ZIV_ID_SEQ geno_oce, ziv.ZIV_OCE_SEQ stari_stars, geno.GEN_CHIP from (select ziv.ZIV_OCE_SEQ from GENOTIPIZIRANE_ZIVALI gen, GOVEDO.ZIVALI ziv where ziv.ziv_id_seq=gen.ZIV_ID_SEQ) oce, GOVEDO.ZIVALI ziv, GENOTIPIZIRANE_ZIVALI geno where oce.ZIV_OCE_SEQ=ziv.ZIV_ID_SEQ and ziv.ZIV_OCE_SEQ=geno.ZIV_ID_SEQ;
 
 select * from govedo.zivali ziv, GENOTIPIZIRANE_ZIVALI gen where gen.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ and gen.ZIV_ID_SEQ=1299471;
 
select distinct potomci.ziv_oce_seq, ziva.ZIV_ID_SEQ from (select ziv.ziv_id_seq, ziv.ziv_oce_seq from govedo.zivali ziv where ziv.ZIV_OCE_SEQ=1299471) potomci, govedo.zivali ziva where potomci.ziv_oce_seq=ziva.ZIV_oce_SEQ ;



create or replace function poisci_gen_stStarse as
begin
select ocetje.zival zival, ocetje.oce oce,ziva.ZIV_ID_SEQ Oce1, ziva.ZIV_OCE_SEQ stari_oce, geno.GEN_CHIP stari_oce_cip from (select distinct ziv.ZIV_ID_SEQ zival, ziv.ZIV_OCE_SEQ oce from govedo.zivali ziv, GENOTIPIZIRANE_ZIVALI gen where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ) ocetje, govedo.zivali ziva, GENOTIPIZIRANE_ZIVALI geno where ocetje.oce=ziva.ZIV_ID_SEQ and ziva.ZIV_OCE_SEQ=geno.ZIV_ID_SEQ
end;

select * from govedo.zivali ziv where ziv.ziv_id_seq=3359948;


select * from govedo.zivali ziv where ziv.STEV_ORIG_ZIVAL=34778803;
