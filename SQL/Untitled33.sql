select * from GOVEDO.BIKOVSKE_MATERE bm where extract(year from bm.DAT_STATUS_BM)>2018;


SELECT
  COUNT (DISTINCT gen.ZIV_ID_SEQ),
  ziv.SIF_SPOL
FROM
  GOVEDO.ZIVALI_PASMA_SPADA pasma,
  GENOTIPIZIRANE_ZIVALI gen,
  zivali ziv
WHERE
  pasma.ZIV_ID_SEQ   =ziv.ZIV_ID_SEQ
AND gen.ZIV_ID_SEQ   =ziv.ZIV_ID_SEQ
AND pasma.PASMA_SPADA=1
GROUP BY
  ziv.SIF_SPOL;
  
  
select * from GENOTIPIZIRANE_ZIVALI gen where gen.ZIV_ID_SEQ=3963934;

select ziv.ZIV_ID_SEQ, stars.MATI_SEQ, stars.OCE_SEQ, ziv.SP1_SIFRA_PASMA, ziv.SIF_SPOL,stars.STEVILO_NEUJEMANJ_KONCNO from GENOTIPIZIRANE_ZIVALI gen, zivali ziv, SNP_PREVERJANJE_STARSEVSTVO stars where gen.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ and stars.ZIVAL_SEQ=ziv.ZIV_ID_SEQ and stars.STEVILO_NEUJEMANJ_KONCNO > 300;


select * from zivali where ZIV_ID_SEQ=4025718;
select * from janao.genotipizirane_zivali gen where  gen.ziv_id_seq=4025718;


select * from zivali where STEV_ORIG_ZIVAL='14819108';

select * from PARENTAL_SNP800 snp where snp.ZGP_ZIV_ID_SEQ=4025718 and not snp.ZGP_VREDNOST = ('A') and not snp.ZGP_VREDNOST = ('B');
select * from PARENTAL_SNP800 snp where snp.ZGP_ZIV_ID_SEQ=1603628 and not snp.ZGP_VREDNOST = ('A') and not snp.ZGP_VREDNOST = ('B');