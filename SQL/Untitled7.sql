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
   FROM THE (SELECT CAST (govedo.izpis_prireja_mleka.zivali_lok_obdobje (
                             20029,
                             to_DATE('01.01.1990', 'dd.mm.yyyy'),
                             TO_DATE('01.01.2014', 'dd.mm.yyyy'),
1) AS govedo.ziv_lok_obd_table)
               FROM DUAL)
--ORDER BY govedo_splosne_procedure.sortiraj_stevilko_zivali (drzava, stevilka)
;

SELECT *
 FROM THE (
         SELECT CAST (
                   govedo.kontrola_mlecnosti_izracun.f_izracun_laktacije_obdobje (
                      20029,
                      TO_DATE ('01.01.12', 'dd.mm.yyyy'),
                             TO_DATE ('31.12.12', 'dd.mm.yyyy')) as govedo.laktacije_table)
           FROM DUAL);