select ziv.DRZ_ORIG_ZIVAL, ziv.STEV_ORIG_ZIVAL, ziv.DAT_ROJSTVO, ziv.SIF_SPOL, gen.ZIV_ID_SEQ, gen.GEN_CHIP, gen.GEN_DATUM from GENOTIPIZIRANE_ZIVALI gen, govedo.zivali ziv where ziv.ziv_id_seq=gen.ZIV_ID_SEQ ;

select * from govedo.zivali ziv where ziv.STEV_ORIG_ZIVAL='04231516';


select * from govedo.zivali ziv where ziv.ZIV_ID_SEQ=347112;

select distinct gen.ZIV_ID_SEQ from GENOTIPIZIRANE_ZIVALI gen, GOVEDO.ZIVALI_PASMA_SPADA pasma, govedo.zivali ziv where ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ and pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ;



truncate table PARENTAL_SNP800_TEMP;



insert into PARENTAL_SNP800 par
(par.ZGP_ZIV_ID_SEQ,par.ZGP_SGPL_SIFRA,par.ZGP_vrednost)
select partemp.ZGP_ZIV_ID_SEQ, partemp.ZGP_SGPL_SIFRA, partemp.ZGP_VREDNOST 
from PARENTAL_SNP800_TEMP partemp;


truncate table GENOTIPIZIRANE_ZIVALI_TEMP;

INSERT INTO GENOTIPIZIRANE_ZIVALI gen
(gen.ZIV_ID_SEQ,gen.GEN_CHIP, gen.GEN_DATUM )
SELECT temp.ZIV_ID_SEQ, temp.GEN_CHIP, temp.GEN_DATUM
FROM GENOTIPIZIRANE_ZIVALI_TEMP temp;


select * from govedo.zivali ziv where ziv.ZIV_ID_SEQ=4516927;

--bikovske matere
select distinct ziv.ZIV_ID_SEQ, ziv.DRZ_ORIG_ZIVAL  || ziv.STEV_ORIG_ZIVAL, mat.SIF_USMER_BM, pasma.PASMA_SPADA from GOVEDO.BIKOVSKE_MATERE mat, govedo.zivali ziv, GOVEDO.ZIVALI_PASMA_SPADA pasma where ziv.ZIV_ID_SEQ=mat.BM_ziv_ID_SEQ and pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ and (pasma.PASMA_SPADA=2 OR pasma.PASMA_SPADA=3 OR pasma.PASMA_SPADA=1);

DELETE FROM PARENTAL_SNP800 snp
WHERE snp.ZGP_ZIV_ID_SEQ IN (SELECT temp.ZGP_ZIV_ID_SEQ FROM PARENTAL_SNP800_TEMP temp );

DELETE FROM GENOTIPIZIRANE_ZIVALI gen
WHERE gen.ZIV_ID_SEQ IN (SELECT temp.ZIV_ID_SEQ FROM GENOTIPIZIRANE_ZIVALI_TEMP temp );



select * from SNP_PREVERJANJE_STARSEVSTVO por where por.ZIVAL_SEQ=4516927;
select * from GENOTIPIZIRANE_ZIVALI por where por.ZIV_ID_SEQ=4516927;
select * from PARENTAL_SNP800 par where par.ZGP_ZIV_ID_SEQ=4516927;


select count(*) from PARENTAL_SNP800 par,
(select distinct gen.ZIV_ID_SEQ zival, gen.GEN_DATUM from GENOTIPIZIRANE_ZIVALI gen, govedo.zivali ziv, GOVEDO.ZIVALI_PASMA_SPADA pasma where gen.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ and ziv.ZIV_ID_SEQ=pasma.ZIV_ID_SEQ and pasma.PASMA_SPADA=1) rjava where rjava.zival=par.ZGP_ZIV_ID_SEQ;


select * from govedo.zivali ziv where ziv.ZIV_ID_SEQ=003933473;

select gen.ZIV_ID_SEQ from GENOTIPIZIRANE_ZIVALI gen
minus
select par.ZGP_ZIV_ID_SEQ from PARENTAL_SNP800 par;


select count(distinct (por.ZIVAL_SEQ)) from SNP_POREKLO_PREVERJANJE por;
select count(distinct (gen.ZIV_ID_SEQ)) from GENOTIPIZIRANE_ZIVALI gen;


select count(*) from SNP_PREVERJANJE_STARSEVSTVO por where por.STEVILO_NEUJEMANJ_KONCNO=0;