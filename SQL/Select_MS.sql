drop table gp_jana

create table gp_jana as
SELECT  DISTINCT DRZ_ORIG_ZIVAL
      
      || STEV_ORIG_ZIVAL ID,
      ZGP_ZIV_ID_SEQ,
      SIF_SPOL SPOL,
      VREDNOST1,
      VREDNOST2,
      SGPL_NAZIV
    FROM ANDREJABO.ANDREJA_IZPIS_MS,
      govedo.zivali
    WHERE ZGP_ZIV_ID_SEQ=ZIV_ID_SEQ
   -- AND ZIV_ID_SEQ in ()                  tukaj navedi sekvence (lahko nadomestiš z ID ji STEV_ORIG_ZIVAL)
AND govedo.zivali.SP1_SIFRA_PASMA=1
and  DRZ_ORIG_ZIVAL='SI'
    order by ZGP_ZIV_ID_SEQ;
 
SELECT  seznam_ziv.ID,
nvl(p1.vrednost1,'?') as BM1824, nvl(p1.vrednost2,'?') as BM1824_1,
nvl(p2.vrednost1,'?') as BM2113 ,nvl(p2.vrednost2,'?') as BM2113_1 ,
nvl(p3.vrednost1,'?') as ETH10, nvl(p3.vrednost2,'?') as ETH10_1,
nvl(p4.vrednost1,'?') as  ETH225, nvl(p4.vrednost2,'?') as  ETH225_1,
nvl(p5.vrednost1,'?') as ETH3, nvl(p5.vrednost2,'?') as ETH3_1,
nvl(p6.vrednost1,'?') as INRA23,nvl(p6.vrednost2,'?') as INRA23_1,
nvl(p7.vrednost1,'?') as SPS115, nvl(p7.vrednost2,'?') as SPS115_1,
nvl(p8.vrednost1,'?') as TGLA122, nvl(p8.vrednost2,'?') as TGLA122_1,
nvl(p9.vrednost1,'?') as TGLA126, nvl(p9.vrednost2,'?') as TGLA126_1,
nvl(p10.vrednost1,'?') as TGLA227, nvl(p10.vrednost2,'?') as TGLA227_1,
nvl(p11.vrednost1,'?') as TGLA53,nvl(p11.vrednost2,'?') as TGLA53_1,
nvl(p12.vrednost1,'?') as BM1818,nvl(p12.vrednost2,'?') as BM1818_1
 
 
  FROM   (select ID, sgpl_naziv,   vrednost1,   vrednost2 from janao.gp_jana where sgpl_naziv='BM1824' ) p1,
         (select ID, sgpl_naziv,  vrednost1,  vrednost2 from janao.gp_jana where sgpl_naziv='BM2113' ) p2,
         (select ID, sgpl_naziv,  vrednost1,  vrednost2 from janao.gp_jana where sgpl_naziv='ETH10' ) p3,
         (select ID, sgpl_naziv,  vrednost1,  vrednost2 from janao.gp_jana where sgpl_naziv='ETH225' ) p4,
         (select ID, sgpl_naziv,  vrednost1,  vrednost2 from janao.gp_jana where sgpl_naziv='ETH3' ) p5,
         (select ID, sgpl_naziv,  vrednost1,  vrednost2 from janao.gp_jana where sgpl_naziv='INRA23' ) p6,
         (select ID, sgpl_naziv,  vrednost1,  vrednost2 from janao.gp_jana where sgpl_naziv='SPS115' ) p7,
         (select ID, sgpl_naziv,  vrednost1,  vrednost2 from janao.gp_jana where sgpl_naziv='TGLA122' ) p8,
         (select ID, sgpl_naziv,  vrednost1,  vrednost2 from janao.gp_jana where sgpl_naziv='TGLA126' ) p9,
         (select ID, sgpl_naziv,  vrednost1,  vrednost2 from janao.gp_jana where sgpl_naziv='TGLA227' ) p10,
         (select ID, sgpl_naziv,  vrednost1,  vrednost2 from janao.gp_jana where sgpl_naziv='TGLA53' ) p11,
         (select ID, sgpl_naziv,  vrednost1,  vrednost2 from janao.gp_jana where sgpl_naziv='BM1818' ) p12,
   
        
        
 (select distinct ID 
 from JANAO.gp_jana) seznam_ziv
WHERE  
seznam_ziv.ID = p1.ID (+) and
seznam_ziv.ID = p2.ID (+) and
seznam_ziv.ID = p3.ID (+) and
seznam_ziv.ID = p4.ID (+) and
seznam_ziv.ID = p5.ID (+) and
seznam_ziv.ID = p6.ID (+) and
seznam_ziv.ID = p7.ID (+) and
seznam_ziv.ID = p8.ID (+) and
seznam_ziv.ID = p9.ID (+) and
seznam_ziv.ID = p10.ID (+) and
seznam_ziv.ID = p11.ID (+) and
seznam_ziv.ID = p12.ID (+)
 
 
order by
 seznam_ziv.ID;