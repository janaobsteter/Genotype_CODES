SELECT zival_seq,
          mati_seq,
          oce_seq,
          REPLACE (lastnost, '_', '') lastnost,
          zival_vrednost1,
          zival_vrednost2
     FROM (SELECT zsgpl_1.zgp_ziv_id_seq zival_seq,
                  zsgpl_1.ziv_mati_seq mati_seq,
                  zsgpl_1.ziv_oce_seq oce_seq,
                  SUBSTR (zsgpl_1.sgpl_naziv,
                          1,
                          INSTR (zsgpl_1.sgpl_naziv, '_'))
                     lastnost,
                  zsgpl_1.zgp_vrednost zival_vrednost1,
                  zsgpl_2.zgp_vrednost zival_vrednost2
             FROM (SELECT zgp_ziv_id_seq,
                          ziv_mati_seq,
                          ziv_oce_seq,
                          sgpl_naziv,
                          zgp_vrednost
                     FROM janao.sifrant_SNP800, janao.PARENTAL_SNP800, JANAO.NAPACNE_MAME
                    WHERE     SUBSTR (sgpl_naziv, -2) = '_1'
                          AND sgpl_sifra = zgp_sgpl_sifra
                          AND zgp_ziv_id_seq = ziv_id_seq) zsgpl_1,
                  (SELECT zgp_ziv_id_seq, sgpl_naziv, zgp_vrednost
                     FROM janao.sifrant_SNP800, janao.PARENTAL_SNP800
                    WHERE SUBSTR (sgpl_naziv, -2) = '_2'
                          AND sgpl_sifra = zgp_sgpl_sifra) zsgpl_2
            WHERE zsgpl_1.zgp_ziv_id_seq = zsgpl_2.zgp_ziv_id_seq
                  AND SUBSTR (zsgpl_1.sgpl_naziv,
                              1,
                              INSTR (zsgpl_1.sgpl_naziv, '_')) =
                         SUBSTR (zsgpl_2.sgpl_naziv,
                                 1,
                                 INSTR (zsgpl_2.sgpl_naziv, '_')));