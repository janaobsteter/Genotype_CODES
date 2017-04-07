 
SELECT  DISTINCT DRZ_ORIG_ZIVAL
      ||' '
      || STEV_ORIG_ZIVAL ID,
      ZGP_ZIV_ID_SEQ,
      SIF_SPOL SPOL,
      VREDNOST1,
      VREDNOST2,
      SGPL_NAZIV
    FROM ANDREJABO.AN,
      govedo.zivali
    WHERE ZGP_ZIV_ID_SEQ=ZIV_ID_SEQ
   -- AND ZIV_ID_SEQ in ()                  tukaj navedi sekvence (lahko nadomestiš z ID ji STEV_ORIG_ZIVAL)
--AND govedo.zivali.SP1_SIFRA_PASMA=5      lahko iščeš samo po določeni pasmi
    order by ZGP_ZIV_ID_SEQ;