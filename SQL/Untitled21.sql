INSERT INTO GENOTIPIZIRANE_ZIVALI
(ziv_id_Seq, GEN_CHIP,  GEN_DATUM)
SELECT ziv_id_Seq, GEN_CHIP, GEN_DATUM
FROM GENOTIPIZIRANE_ZIVALI_TEMP

INSERT INTO janao.PARENTAL_SNP800
(ZGP_ziv_id_Seq, ZGP_SGPL_SIFRA, ZGP_VREDNOST)
SELECT ZGP_ziv_id_Seq, ZGP_SGPL_SIFRA, ZGP_VREDNOST
FROM janao.parental_SNP800_TEMP;

select ziv.ZIV_ID_SEQ, ziv.DRZ_ORIG_ZIVAL, ziv.STEV_ORIG_ZIVAL, ziv.DAT_ROJSTVO from govedo.zivali ziv where ziv.SP1_SIFRA_PASMA=1;

select distinct par.ZGP_ZIV_ID_SEQ from PARENTAL_SNP800_TEMP par;

select distinct par.ZIVAL_SEQ, ZIV.SP1_SIFRA_PASMA from SNP_PREVERJANJE_STARSEVSTVO par, GOVEDO.ZIVALI ZIV WHERE ZIV.ZIV_ID_SEQ=par.ZIVAL_SEQ;


 alter table   GENOTIPIZIRANE_ZIVALI drop column ID_ZIVALI ;
   
UPDATE PARENTAL_SNP800 par SET par.ZGP_VREDNOST = REPLACE(par.ZGP_VREDNOST, 'NA', '');


select distinct par.ZGP_VREDNOST from PARENTAL_SNP800_TEMP par;

select distinct ZIV_ID_SEQ from GENOTIPIZIRANE_ZIVALI;

select distinct ZIVAL_SEQ from SNP_PREVERJANJE_STARSEVSTVO;




select * from GENOTIPIZIRANE_ZIVALI gen where gen.GEN_CHIP="GGPv04";


select ms.ZIVAL_SEQ, sum(ms.PREVERJANJE_KONCNO) from GOVEDO.POREKLO_PREVERJANJE ms, SNP_POREKLO_PREVERJANJE snp where ms.ZIVAL_SEQ=snp.ZIVAL_SEQ and ms.PREVERJANJE_KONCNO is not null group by ms.ZIVAL_SEQ;

select distinct * from GOVEDO.GEN_RJAVA gen, GENOTIPIZIRANE_ZIVALI jana where jana.ZIV_ID_SEQ=gen.GEN_ZIV_ID_SEQ;


select ziv.ZIV_ID_SEQ, ziv.DRZ_ORIG_ZIVAL, ziv.STEV_ORIG_ZIVAL, ziv.DRZ_ORIG_ZIVAL || ziv.STEV_ORIG_ZIVAL  from GOVEDO.ZIVALI ziv where  ziv.SP1_SIFRA_PASMA = 1;


--najdi brate genotipiziranih bikov, ki so tudi genotipizirani /ali pa očete
select ziv.ZIV_ID_SEQ, ziv.ZIV_OCE_SEQ, gen.GEN_CHIP, gen.ZIV_ID_SEQ genID from GENOTIPIZIRANE_ZIVALI gen, GOVEDO.ZIVALI ziv where ziv.ZIV_OCE_SEQ=gen.ZIV_ID_SEQ and ziv.ZIV_ID_SEQ=gen.ZIV_ID_SEQ;



select ziv.ZIV_ID_SEQ, ziv.ZIV_MATI_SEQ, ziv.ZIV_OCE_SEQ from govedo.zivali  ziv where ziv.ZIV_ID_SEQ=4077913;
