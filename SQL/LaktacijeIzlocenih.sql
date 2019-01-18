select * from (
SELECT
  izlocene.seqZivali,
  izlocene.pasma,
  izlocene.rojstvo,
  izlocene.PGK_DAT_ZAKOLA,
  izlocene.PGK_TOPLA_MASA,
  izlocene.PGK_STAROST_DNI,
  izlocene.telitev_izlocitve,
  teli.DAT_TELITEV,
  izlocene.pasma_spada,
  lak.KG_MLEKO_305
FROM
  (
    SELECT DISTINCT
      ziv.ziv_id_seq seqZivali,
      ziv.SP1_SIFRA_PASMA pasma,
      ziv.DAT_ROJSTVO rojstvo,
      klav.PGK_DAT_ZAKOLA,
      klav.PGK_TOPLA_MASA,
      klav.PGK_STAROST_DNI,
      pasmaSpada.PASMA_SPADA,
      MAX(tel.ZAP_TELITEV) telitev_izlocitve
    FROM
      govedo.zivali ziv,
      GOVEDO.POD_GOV_KLAVNICE klav,
      GOVEDO.TELITVE tel,
      GOVEDO.ZIVALI_PASMA_SPADA pasmaSpada
    WHERE
      ziv.ZIV_ID_SEQ  =klav.PGK_ZIV_ID_SEq
    AND ziv.ZIV_ID_SEQ=tel.TEL_ZIV_ID_SEQ
    and pasmaSpada.ZIV_ID_SEQ=ziv.ZIV_ID_SEQ
    GROUP BY
      ziv.ziv_id_seq,
      klav.PGK_DAT_ZAKOLA,
      ziv.DAT_ROJSTVO,
      klav.PGK_TOPLA_MASA,
      pasmaSpada.PASMA_SPADA,
      klav.PGK_STAROST_DNI,
      ziv.SP1_SIFRA_PASMA
  )
  izlocene, --tabela zadnjih laktacij izločenih živali
  GOVEDO.TELITVE teli,
  GOVEDO.LAKTACIJE lak
WHERE
  izlocene.seqZivali          =teli.TEL_ZIV_ID_SEQ
  and teli.TEL_ID_SEQ=lak.LAK_TEL_ID_SEQ (+)
AND izlocene.telitev_izlocitve=teli.ZAP_TELITEV) zadnje_lak where zadnje_lak.KG_MLEKO_305 is not null;