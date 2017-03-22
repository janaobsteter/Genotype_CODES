select *
from GOVEDO.LOKACIJE lok
where lok.KMG_MID=100253542 --pridobi številko črede
;


--zberi vse živali, tudi tiste s selitvami
--osnovni podatki
--poglej, koliko je krav in koliko telic
--ne samo telitev na tej lokaciji
--drop table vse_zivali_Gasper_Franc
create table vse_zivali_Gasper_Franc  as (
select 
  ziv.ZIV_ID_SEQ,
  ziv.DRZ_ORIG_ZIVAL
  || ziv.STEV_ORIG_ZIVAL AS ID_zivali,
  pasma_sp.PASMA_SPADA as pasma_spada,
  pasma.IME_PASMA_SLO,
  ziv.IME_ZIVAL,
  ziv.DAT_ROJSTVO,
  ziv.SIF_SPOL,
  max(tel.ZAP_TELITEV) tel
from GOVEDO.ZIVALI ziv,
  GOVEDO.ZIVALI_PASMA_SPADA pasma_sp,
  GOVEDO.SIFRANT_PASEM pasma,
  GOVEDO.TELITVE tel,
  (select ziv.ZIV_ID_SEQ as ziv_ID
  from GOVEDO.ZIVALI ziv
  where ziv.CRE_SIFRA_CREDA=20029
    UNION --union=distinct results, union all = all rows
  select selitve.OPZ_ZIV_ID_SEQ as ziv_ID
  from GOVEDO.ODJAVA_PRIJAVA_ZIVALI selitve
  where selitve.CRE_NOVA_SIFRA_CREDA=20029) tab1
where ziv.ZIV_ID_SEQ=tab1.ziv_ID
and pasma_sp.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
and pasma_sp.PASMA_SPADA=pasma.SIFRA_PASMA(+)
and ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ(+)
group by ziv.ZIV_ID_SEQ,
  ziv.DRZ_ORIG_ZIVAL,
  ziv.STEV_ORIG_ZIVAL,
  pasma_sp.PASMA_SPADA,
  pasma.IME_PASMA_SLO,
  ziv.IME_ZIVAL,
  ziv.DAT_ROJSTVO,
  ziv.SIF_SPOL
);
--1553 ŽIVALI SKUPNO


--preveri, koliko telitev ima NULL lokacijo telitve IN lokacijo presušitve IN lokacijo kontrole --> 0
select distinct vse_zivali.ZIV_ID_SEQ--, tel.TEL_LOKACIJA_TELITVE, lak.LOK_PRESUSITEV, kontrole.SIFRA_LOKACIJE
from VSE_ZIVALI_GASPER_FRANC vse_zivali, GOVEDO.TELITVE tel, GOVEDO.LAKTACIJE lak, GOVEDO.KONTROLA_MLECNOSTI kontrole
where vse_zivali.ZIV_ID_SEQ= tel.TEL_ZIV_ID_SEQ and tel.TEL_ID_SEQ=lak.LAK_TEL_ID_SEQ and kontrole.KM_TEL_ID_SEQ=tel.TEL_ID_SEQ
and tel.TEL_LOKACIJA_TELITVE is null --and lak.LOK_PRESUSITEV is null --and kontrole.SIFRA_LOKACIJE is null;



--select vse živali z lokacijo telitve / lokacijo presušitve / kontrolo laktacije na lokaciji 20029
--drop table Gasper_Franc_telitve
create table Gasper_Franc_telitve as (
select distinct vse_zivali.ID_ZIVALI, vse_zivali.ZIV_ID_SEQ, tel.tel_ID_SEQ, tel.TEL_LOKACIJA_TELITVE, lak.LOK_PRESUSITEV
from VSE_ZIVALI_GASPER_FRANC vse_zivali, GOVEDO.KONTROLA_MLECNOSTI kontrola, GOVEDO.LAKTACIJE lak, GOVEDO.TELITVE tel
where (tel.TEL_LOKACIJA_TELITVE=20029  OR lak.LOK_PRESUSITEV=20029 OR kontrola.SIFRA_LOKACIJE=20029 )
 and tel.TEL_ZIV_ID_SEQ=vse_zivali.ZIV_ID_SEQ
 and lak.LAK_TEL_ID_SEQ=tel.TEL_ID_SEQ 
 and kontrola.KM_TEL_ID_SEQ=tel.TEL_ID_SEQ
);

--tabela vseh laktacij in telitev (data v R) 
--podatki o laktacijah - maksimalna st. cela, maks količina belja in maščob na leto
--tukaj črpaš seznam živali (sekvenco) iz kreirane tabele zivali_telitve (vključiš samo tiste, ki so telile na tej lokaciji)
SELECT 
  ziv.ZIV_ID_SEQ,
  ziv.DRZ_ORIG_ZIVAL || ziv.STEV_ORIG_ZIVAL as ID_zivali,
  tel.TEL_ID_SEQ,
  tel.ZAP_TELITEV,
  tel.DAT_TELITEV,
  lak.KG_MLEKO_305,
  lak.KG_MAST_305,
  lak.KG_BELJAK_305,
  lak.KG_MLEKO_CELA_LAKT
FROM 
  GASPER_FRANC_TELITVE GF_tel,
  GOVEDO.ZIVALI ziv,
  GOVEDO.LAKTACIJE lak,
  GOVEDO.TELITVE tel 
WHERE 
ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ
and tel.TEL_ID_SEQ=GF_tel.TEL_ID_SEQ
and tel.TEL_ID_SEQ=lak.LAK_TEL_ID_SEQ(+)
;



--ZDRUŽI PO ŽIVALIH S TELITVAMI
--tabela vseh distinct živali s telitvami (vsa v R)
--podatki o vsoti celih laktacij (življenjska prireja), povprečni celi in standardni laktaciji živali, ki so identificirano (GF_tel) imele laktacije na tej lokaciji
--povprečna količina maščob, beljakovin in maščob + beljakovin v st. laktaciji
SELECT 
  ziv.ZIV_ID_SEQ,
  ziv.DRZ_ORIG_ZIVAL
  || ziv.STEV_ORIG_ZIVAL AS ID_zivali,
  pasma_sp.PASMA_SPADA as pasma_spada,
  pasma.IME_PASMA_SLO,
  ziv.IME_ZIVAL,
  ziv.DAT_ROJSTVO,
  ziv.SIF_SPOL,
  (SELECT pv.VREDNOST_12_PV FROM ARHIV.PLEMENSKE_VREDNOSTI pv WHERE pv.SIFRA_LAST=79 AND pv.PV_ZIV_ID_SEQ(+)=ziv.ZIV_ID_SEQ) AS SSI,
  round(SUM(lak.KG_MLEKO_CELA_LAKT),2) as MlekoZiv,
  round(avg(lak.KG_MLEKO_CELA_LAKT),2) as lakt_cela_KG_povprecje,
  round(avg(lak.KG_MLEKO_305),2) as mleko_305_KG_povprecje,
  round(avg(lak.KG_MAST_305),2) as mast_305_KG_povprecje,
  round(avg(lak.KG_MAST_305)*100/avg(lak.KG_MLEKO_305), 2) as odstotek_masc_povp,
  round(avg(lak.KG_BELJAK_305),2) as beljak_305_KG_povprecje,
  round(avg(lak.KG_BELJAK_305)*100/avg(lak.KG_MLEKO_305), 2) as odstotek_belj_povp,
  round(avg(lak.KG_MAST_305+lak.KG_BELJAK_305), 2) mast_belj_305__KG_povprecje,
  max (tel.ZAP_TELITEV) st_telitev
FROM GOVEDO.ZIVALI ziv,
  GOVEDO.TELITVE tel,
  GOVEDO.LAKTACIJE lak,
  GOVEDO.ZIVALI_PASMA_SPADA pasma_sp,
  GOVEDO.SIFRANT_PASEM pasma,
  GASPER_FRANC_TELITVE GF_tel
WHERE gf_tel.TEL_ID_SEQ=tel.TEL_ID_SEQ
and ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ
AND tel.TEL_ID_SEQ       =lak.LAK_TEL_ID_SEQ(+)
and pasma_sp.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
and pasma_sp.PASMA_SPADA=pasma.SIFRA_PASMA(+) --poveži tabelo pasma spada z imeni - vendar 222 ni v tabeli imen (+)
GROUP BY ZIV.ZIV_ID_SEQ,
  ziv.DRZ_ORIG_ZIVAL,
  ziv.STEV_ORIG_ZIVAL,
  pasma_sp.PASMA_SPADA,
  pasma.IME_PASMA_SLO,
  ziv.IME_ZIVAL,
  ziv.DAT_ROJSTVO,
  ziv.SIF_SPOL;


--poglej, če je kakšen plemenski bik
select pb.PB_ZIV_ID_SEQ, pasma_sp.PASMA_SPADA, pasma.IME_PASMA_SLO
from GOVEDO.PLEMENSKI_BIKI pb,
  govedo.ZIVALI_PASMA_SPADA pasma_sp,
 GOVEDO.SIFRANT_PASEM pasma,
(select ziv.ZIV_ID_SEQ as ziv_ID
  from GOVEDO.ZIVALI ziv
  where ziv.CRE_SIFRA_CREDA=20029
    UNION --union=distinct results, union all = all rows
  select selitve.OPZ_ZIV_ID_SEQ as ziv_ID
  from GOVEDO.ODJAVA_PRIJAVA_ZIVALI selitve
  where selitve.CRE_NOVA_SIFRA_CREDA=20029) tab1
where tab1.ziv_ID=pb.PB_ZIV_ID_SEQ
and pb.PB_ZIV_ID_SEQ=pasma_sp.ZIV_ID_SEQ
and pasma_sp.PASMA_SPADA=pasma.SIFRA_PASMA
;






