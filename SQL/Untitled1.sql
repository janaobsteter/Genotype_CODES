
       SELECT ziv_id_seq,
         drzava,
         stevilka,
         ime,
         spol,
         rojstvo,
         id_kmetije,
         od_dne,
         do_dne,
         (do_dne - od_dne) kd
    FROM THE (SELECT CAST (izpis_prireja_mleka.zivali_lok_obdobje (
                                    -- 20029
                              &v_id_lok
                              ,
                             -- TO_DATE ('01.01.12', 'dd.mm.rr') 
                              &v_od_dne
                              ,
                             -- TO_DATE ('31.12.12', 'dd.mm.rr')
                              &v_do_dne
                              ,
                              --0
                              &v_vse_zivali
                              )
                              AS ziv_lok_obd_table)
                FROM DUAL)
ORDER BY govedo_splosne_procedure.sortiraj_stevilko_zivali (drzava, stevilka)
;
