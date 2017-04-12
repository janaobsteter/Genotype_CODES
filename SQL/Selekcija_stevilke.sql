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

--število moških in ženskih telet rojenih op letih, povprečje 2010-2016
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
    AND extract(YEAR FROM tel.DAT_TELITEV) BETWEEN 2010 AND 2016
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
    AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2010 AND 2016
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
    AND extract(YEAR FROM ziv.DAT_ROJSTVO) BETWEEN 2010 AND 2016
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
    AND extract(YEAR FROM ziv.DAT_ROJSTVO)  in (2014,2015,2016) -- between 2000 AND 2016
    and ziv.SIF_SPOL=1
    and ziv.ZIV_ID_SEQ=pb.PB_ZIV_ID_SEQ
    and ss.SS_ID_SEQ=17
    and pasma.PASMA_SPADA=1
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
      
      
--genotipizirane rojene 2015      
select count(distinct gen.ZIV_ID_SEQ) from JANAO.ALL_GEN_IND_27012017 gen, govedo.zivali ziv, GOVEDO.ZIVALI_PASMA_SPADA PASMA
where gen.ZIV_ID_SEQ=ziv.ziv_id_seq
AND ziv.ZIV_ID_SEQ=PASMA.ZIV_ID_SEQ
AND ziv.SIF_SPOL=1
AND PASMA.PASMA_SPADA=1
and (extract(year from ziv.DAT_ROJSTVO))=2015;

--koliko potomk po bikih v pripustu, koliko po AI
select  count(distinct ziv.ZIV_ID_SEQ), extract(YEAR FROM ziv.DAT_ROJSTVO) DatRoj, pb.SIF_UPORABA_PB from  govedo.zivali ziv, GOVEDO.ZIVALI_PASMA_SPADA pasma , GOVEDO.PLEMENSKI_BIKI pb,
      GOVEDO.SEZNAM_SIFRANTOV ss,
      GOVEDO.ZBIRKA_SIFRANTOV_GOVEDO zs
where ziv.ziv_id_seq=pasma.ZIV_ID_SEQ and pasma.PASMA_SPADA=1 --and ziv.SIF_SPOL=2
and       zs.ZSG_SS_ID_SEQ=ss.SS_ID_SEQ and 
      pb.SIF_UPORABA_PB=zs.ZSG_SIFRA
    AND extract(YEAR FROM ziv.DAT_ROJSTVO) in(2014) --BETWEEN 2000 AND 2016 -
    and ziv.ZIV_OCE_SEQ=pb.PB_ZIV_ID_SEQ
    and ss.SS_ID_SEQ=17
    group by pb.SIF_UPORABA_PB,extract(YEAR FROM ziv.DAT_ROJSTVO) ;
    
    
--potomke rojene v dolocenem letu --> koliko ocetov iz pripsta, AI, genom. -_> st. ocetov
select  count(distinct ziv.ZIV_OCE_SEQ) st_potomcev,  pb.SIF_UPORABA_PB from  govedo.zivali ziv, GOVEDO.ZIVALI_PASMA_SPADA pasma , GOVEDO.PLEMENSKI_BIKI pb
where ziv.ziv_id_seq=pasma.ZIV_ID_SEQ and pasma.PASMA_SPADA=1 --and ziv.SIF_SPOL=2
    and ziv.ZIV_OCE_SEQ=pb.PB_ZIV_ID_SEQ
    AND (extract(YEAR FROM ziv.DAT_ROJSTVO)) in(2014) --BETWEEN 2000 AND 2016 -
    and count(distinct ziv.ZIV_OCE_SEQ) >= 50
    group by pb.SIF_UPORABA_PB;


--koliko očetov v AI več kot 50 potomk v letu --> tukaj dobi samo seznam očetov, potem v Ru dobi ven več kot 50 in povprečje za tri leta
select  count(distinct ziv.ZIV_OCE_SEQ) st_potomcev, ziv.ZIV_OCE_SEQ, extract(YEAR FROM ziv.DAT_ROJSTVO), pb.SIF_UPORABA_PB  from  govedo.zivali ziv, GOVEDO.ZIVALI_PASMA_SPADA pasma , GOVEDO.PLEMENSKI_BIKI pb
where ziv.ziv_id_seq=pasma.ZIV_ID_SEQ and pasma.PASMA_SPADA=1 --and ziv.SIF_SPOL=2
    and ziv.ZIV_OCE_SEQ=pb.PB_ZIV_ID_SEQ
    AND (extract(YEAR FROM ziv.DAT_ROJSTVO)) in(2014) --BETWEEN 2000 AND 2016 -
group by ziv.ZIV_OCE_SEQ,extract(YEAR FROM ziv.DAT_ROJSTVO), pb.SIF_UPORABA_PB;
    
    
    --ocetje v ai po letih --> pridobi očete, potem pa za te očete pridobi celotno število potomcev --> da dobiš povprečno število doz / oz. potomcev v pripusti 
    --tudi leto rojstva, da vidiš, kako dolgo so v uporabi
 select oce.ocetje,count(distinct ziv.ZIV_ID_SEQ), pb.SIF_UPORABA_PB from 
 (   select ziv.ZIV_OCE_SEQ  ocetje
from zivali ziv
where extract(year from ziv.DAT_ROJSTVO) between 2014 and 2016
and ziv.SP1_SIFRA_PASMA = 1
and ziv.ZIV_OCE_SEQ is not null
) oce, govedo.zivali ziv, GOVEDO.PLEMENSKI_BIKI pb where oce.ocetje=ziv.ZIV_OCE_SEQ and ziv.ZIV_OCE_SEQ=pb.PB_ZIV_ID_SEQ group by pb.SIF_UPORABA_PB, oce.ocetje;


    --koliko let ima potomce
 select distinct oce.ocetje, ziv.ZIV_ID_SEQ, ziv.DAT_ROJSTVO, pb.SIF_UPORABA_PB from 
 (   select ziv.ZIV_OCE_SEQ  ocetje
from zivali ziv
where extract(year from ziv.DAT_ROJSTVO) between 2014 and 2016
and ziv.SP1_SIFRA_PASMA = 1
and ziv.ZIV_OCE_SEQ is not null
) oce, govedo.zivali ziv, GOVEDO.PLEMENSKI_BIKI pb where oce.ocetje=ziv.ZIV_OCE_SEQ and ziv.ZIV_OCE_SEQ=pb.PB_ZIV_ID_SEQ and oce.ocetje=420978;

select * from govedo.zivali ziv where ziv.ZIV_ID_SEQ=420978 ;





--vsi, ki so v tabli plemenski biki, so bili enkrat mladi - razen tisti, ki se začnejo z 8 (pripust) in 7 (uvoženo) - REPUBLISKA
--republiska številka označena kot rodovniška v tabelah
--sicer je sedaj rodovniška številka enaka identifikacijski številki (t.j. SI XXXXXXX - brez SI), včasih pa je bila rodovniška = republiška
select count(distinct plemB.PB_ZIV_ID_SEQ),  (extract (year from ziv.dat_ROJSTVO))
from GOVEDO.PLEMENSKI_BIKI plemB, govedo.zivali ziv
where ziv.ZIV_ID_SEQ=plemB.PB_ZIV_ID_SEQ and ziv.SP1_SIFRA_PASMA=1 
and (extract (year from ziv.dat_ROJSTVO)) between 2006 and 2016
and plemB.STEV_RODOV_PB not like '7%'
and plemB.STEV_RODOV_PB not like '8%'
group by (extract( year from ziv.dat_ROJSTVO))
;

--koliko potomcev rojenih od mladih bikov - to je tri leta po rojstvu bika
select distinct ziv.ZIV_ID_SEQ, oce_pb.ID_PB, oce_pb.ROJ_PB from govedo.zivali ziv, (select distinct plemB.PB_ZIV_ID_SEQ ID_PB,  (extract (year from ziv.dat_ROJSTVO)) ROJ_PB
from GOVEDO.PLEMENSKI_BIKI plemB, govedo.zivali ziv
where ziv.ZIV_ID_SEQ=plemB.PB_ZIV_ID_SEQ and ziv.SP1_SIFRA_PASMA=1 
and (extract (year from ziv.dat_ROJSTVO)) between 2006 and 2016
and plemB.STEV_RODOV_PB not like '7%'
and plemB.STEV_RODOV_PB not like '8%'
) oce_pb
where oce_pb.ID_PB=ziv.ZIV_OCE_SEQ
--and (extract(year from ziv.DAT_ROJSTVO)) between (oce_pb.ROJ_PB) and (oce_pb.ROJ_PB+3)
--group by oce_pb.ID_PB, oce_pb.ROJ_PB
;




--število potomcev po biku šifri 2015
select count(distinct ziv.ziv_id_seq),  pb.SIF_UPORABA_PB from zivali ziv, GOVEDO.PLEMENSKI_BIKI pb  
where (extract(year from ziv.DAT_ROJSTVO))=2015 
and ziv.ZIV_OCE_SEQ=pb.PB_ZIV_ID_SEQ 
and ziv.SP1_SIFRA_PASMA=1
group by pb.SIF_UPORABA_PB;

--koliko je vseh potomcev
select avg(count(distinct ziv.ziv_id_seq)) from zivali ziv where (extract(year from ziv.dat_rojstvo)) between 2014 and 2016 and ziv.SP1_SIFRA_PASMA=1 
group by  (extract(year from ziv.dat_rojstvo));

--še po spolu
select count(distinct ziv.ziv_id_seq), ziv.SIF_SPOL, pb.SIF_UPORABA_PB from zivali ziv, GOVEDO.PLEMENSKI_BIKI pb  
where (extract(year from ziv.DAT_ROJSTVO))=2015 
and ziv.ZIV_OCE_SEQ=pb.PB_ZIV_ID_SEQ 
and ziv.SP1_SIFRA_PASMA=1
group by pb.SIF_UPORABA_PB, ziv.SIF_SPOL;

--koliko je trenutno krav
select count( distinct ziv.ZIV_ID_SEQ) from zivali ziv, GOVEDO.LAKTACIJE lak, GOVEDO.TELITVE tel, GOVEDO.LOKACIJE_STEVILO_ZIVALI lok
where ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ
and tel.TEL_ID_SEQ=lak.LAK_TEL_ID_SEQ
and ziv.AKTIVNA=1
and ziv.SP1_SIFRA_PASMA=1
and ziv.CRE_SIFRA_CREDA=lok.LOKACIJA
and lok.VRSTA_KONTROLE='AP';


--koliko ženskih telic pride do krav
select count(distinct ziv.ZIV_ID_SEQ),(extract(year from ziv.dat_rojstvo)) from zivali ziv, GOVEDO.LOKACIJE_STEVILO_ZIVALI lok, GOVEDO.TELITVE tel
where
ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ
--and tel.TEL_ID_SEQ=lak.LAK_TEL_ID_SEQ
and ziv.SIF_SPOL=2
and ziv.SP1_SIFRA_PASMA=1
and (extract(year from ziv.dat_rojstvo)) between 2010 and 2016
and ziv.CRE_SIFRA_CREDA=lok.LOKACIJA
and lok.VRSTA_KONTROLE='AP'
group by (extract(year from ziv.dat_rojstvo));

--koliko ženskih telic pride do krav
select ziv.ZIV_ID_SEQ,(extract(year from ziv.dat_rojstvo)) from zivali ziv, GOVEDO.LAKTACIJE lak, GOVEDO.TELITVE tel, GOVEDO.LOKACIJE_STEVILO_ZIVALI lok
where ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ
and tel.TEL_ID_SEQ=lak.LAK_TEL_ID_SEQ
and ziv.SP1_SIFRA_PASMA=1
and ziv.CRE_SIFRA_CREDA=lok.LOKACIJA
and lok.VRSTA_KONTROLE='AP';

--koliko je gospodarskih križanj
--koliko je vseh potomcev rjavih krav
select count(distinct potomci.zivID ), POTOMCI.datR
from zivali zivM, 
(select distinct ziv.ziv_id_seq zivID, ziv.ZIV_MATI_SEQ,  (extract(year from ziv.dat_rojstvo)) datR
from zivali ziv 
where (extract(year from ziv.dat_rojstvo)) between 2010 and 2016 
--and ziv.SP1_SIFRA_PASMA=1
) potomci
where potomci.ZIV_MATI_SEQ=zivM.ZIV_ID_SEQ
and zivM.SP1_SIFRA_PASMA=1
group by  POTOMCI.datR;

select count(distinct ziv.ZIV_ID_SEQ), ziv.SIF_SPOL from zivali ziv where ziv.SP1_SIFRA_PASMA=1 and ziv.AKTIVNA=1 group by ziv.SIF_SPOL;
