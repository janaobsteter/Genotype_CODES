
select * from  
(select ziv.ZIV_ID_SEQ biki_id, ziv.DAT_ROJSTVO bikRojstvo
from JANAO.ALL_GEN_IND_27012017 jziv,zivali ziv,rejci_mlekarna
where ziv.ziv_id_seq =jziv.ziv_id_seq
and RM_SIF_LOKACIJA= CRE_SIFRA_CREDA
AND RM_REJEC_ODDAJA_MLEKO = 1
and ziv.SIF_SPOL=1) biki,
govedo.zivali ziv
where biki.biki_id=ziv.ZIV_OCE_SEQ
and ziv.SIF_SPOL=2
and ziv.AKTIVNA=1;


select count(ziv.ZIV_ID_SEQ),extract(year from ziv.DAT_ROJSTVO)
from JANAO.ALL_GEN_IND_27012017 jziv,zivali ziv,rejci_mlekarna
where ziv.ziv_id_seq =jziv.ziv_id_seq
and RM_SIF_LOKACIJA= CRE_SIFRA_CREDA
AND RM_REJEC_ODDAJA_MLEKO = 1
and ziv.AKTIVNA=1
and ziv.SIF_SPOL=2
group by extract(year from ziv.DAT_ROJSTVO)
;

select * from ALL_GEN_IND_27012017 gen, 
(select ziv.ZIV_ID_SEQ seqZiv, ziv.ZIV_OCE_SEQ seqOce, ziv.DAT_ROJSTVO, ziv.AKTIVNA, ziv.SP1_SIFRA_PASMA, ziv.SIF_SPOL from govedo.zivali ziv,GOVEDO.REJCI_MLEKARNA rm
where RM_SIF_LOKACIJA= CRE_SIFRA_CREDA
AND RM_REJEC_ODDAJA_MLEKO = 1) mc_zivali
where gen.ZIV_ID_SEQ=mc_zivali.seqOce;

select distinct mc_zivali.seqOce OCE from ALL_GEN_IND_27012017 gen, 
(select ziv.ZIV_ID_SEQ seqZiv, ziv.ZIV_OCE_SEQ seqOce, ziv.DAT_ROJSTVO, ziv.AKTIVNA, ziv.SP1_SIFRA_PASMA, ziv.SIF_SPOL from govedo.zivali ziv,GOVEDO.REJCI_MLEKARNA rm
where RM_SIF_LOKACIJA= CRE_SIFRA_CREDA
AND RM_REJEC_ODDAJA_MLEKO = 1 and ziv.SIF_SPOL=2 and ziv.AKTIVNA=1) mc_zivali
where gen.ZIV_ID_SEQ=mc_zivali.seqOce;


--očetje teh živali, podatki o pasmi očeta
select pasma.PASMA_SPADA pasmaOce, ocetje.* from GOVEDO.ZIVALI_PASMA_SPADA pasma, 
(select * from ALL_GEN_IND_27012017 gen, 
(select ziv.ZIV_ID_SEQ seqZiv, ziv.ZIV_OCE_SEQ seqOce, ziv.DAT_ROJSTVO, ziv.AKTIVNA, ziv.SP1_SIFRA_PASMA, ziv.SIF_SPOL from govedo.zivali ziv,GOVEDO.REJCI_MLEKARNA rm
where RM_SIF_LOKACIJA= CRE_SIFRA_CREDA
AND RM_REJEC_ODDAJA_MLEKO = 1) mc_zivali
where gen.ZIV_ID_SEQ=mc_zivali.seqOce) ocetje
where ocetje.seqOce=pasma.ZIV_ID_SEQ;

--aktivne ženske
select * from ALL_GEN_IND_27012017 gen, 
(select ziv.ZIV_ID_SEQ seqZiv, ziv.ZIV_OCE_SEQ seqOce, ziv.DAT_ROJSTVO, ziv.AKTIVNA, ziv.SP1_SIFRA_PASMA, ziv.SIF_SPOL from govedo.zivali ziv,GOVEDO.REJCI_MLEKARNA rm
where RM_SIF_LOKACIJA= CRE_SIFRA_CREDA
AND RM_REJEC_ODDAJA_MLEKO = 1 and ziv.SIF_SPOL=2 and ziv.AKTIVNA=1) mc_zivali
where gen.ZIV_ID_SEQ=mc_zivali.seqOce;

SELECT
  pasma.PASMA_SPADA,
  COUNT(*)
FROM
  OCETJEMC,
  govedo.zivali ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma
WHERE
  pasma.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
AND ziv.ZIV_ID_SEQ=OCETJEMC.OCE
GROUP BY
  pasma.pasma_spada;
  
  
  select * from zivali ziv where ziv.ZIV_ID_SEQ in (422733,422148,881061,422336);