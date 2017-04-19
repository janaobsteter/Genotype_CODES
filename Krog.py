# -*- coding: utf-8 -*-
ped = pedigree("~/Documents/PhD/Simulaton/Pedigrees/PedPython.txt")

###################################
# loop za vajo
#TO JE, ČE ŠTARTAŠ IZ NULE - TOREJ POČASI POLNIŠ POPULACIJO!
#################################

stevilo_krogov = 20

for krog in (range(0, stevilo_krogov)):
    if krog == 0:
        ped.set_cat_gen(max(ped.gen), "nr")  # to je samo na prvem loopu
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 == 0], "F")
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 != 0], "M")
        ped.izberi_poEBV_top_catCurrent("F", int(potomciNPn), 'nr', 'potomciNP')
        ped.izberi_poEBV_top_catCurrent("M", int(potomciNPn), 'nr', 'potomciNP')
        categories = ped.save_cat()
        ped = pedigree("~/Documents/PhD/Simulaton/Pedigrees/PedPython.txt")
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 == 0], "F")
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 != 0], "M")
        # prva odbira
        ped.compute_age()
        select_age_0_1(ped)
        ped.add_new_gen_naive(stNB, potomciNPn*2)

        categories.clear()
        categories = ped.save_cat()
        sex = ped.save_sex()

    if krog == 1:
        # SETSEX!!!
        # druga odbira
        ped.set_cat_gen(1, "")
        ped.set_cat_gen(2, "")
        ped.set_cat_old('izl', 'izl', categories)

        ped.compute_age()
        select_age_0_1(ped)
        select_age_1_2(ped)

        ped.add_new_gen_naive(stNB, potomciNPn*2)

        categories.clear()
        categories = ped.save_cat()
        sex = ped.save_sex()
        active = ped.save_active()
        age = ped.save_age()

    if krog >= 2:
        for i in ped.gens():
            ped.set_cat_gen(i, "")

        ped.set_cat_old('izl', 'izl', categories)

        ped.compute_age()
        select_age_0_1(ped)
        select_age_1_2(ped)

        select_age_2_3(ped)

        ped.add_new_gen_naive(stNB, potomciNPn*2)
        ped.compute_age()
        doloci_matere(ped)
        doloci_ocete(ped)
        # ped.set_cat_mother_catCurrent('bm', 'potomciNP') #TO DAJ V FUNKCIJO!

        categories.clear()
        categories = ped.save_cat()
        sex = ped.save_sex()
        active = ped.save_active()
        age = ped.save_age()
        ped.compute_age()  # drugače imaš negativne vrednosti!


#######################################################
#TO JE, ČE ŠTARTAŠ S POLNO AKTIVNO POPULACIJO IN DOLOČIŠ KATEGORIJE
#######################################################
def nastave_cat(PedFile):
    ped = pedigree(PedFile)
    ped.compute_age()
    
    
    #MALES FIRST
    #age 0
    #določi vhlevljene
    ped.izberi_poEBV_top_age_naive(0, vhlevljenin, 'vhlevljeni')
    #določi moška teleta pod 12
    ped.izberi_random_age_naive(0, telMn, 'telM')
    
    #age1
    #določi mlade
    ped.izberi_poEBV_top_age_naive(1, mladin, 'mladi')
    #določi pripust - 1. leto
    ped.izberi_poEBV_OdDo_age_naive(1, mladin, vhlevljenin, 'pripust1')
    #določi bike nad 12 m
    ped.izberi_random_age_naive(1, bik12n, 'bik12')
    
    
    #age2
    ped.izberi_poEBV_top_age_naive(2, mladin, 'cak')
    ped.izberi_poEBV_OdDo_age_naive(1, mladin, (mladin + pripust2n), 'pripust2')
    
    #age3,4
    for i in [3,4]:
        ped.izberi_poEBV_top_age_naive(i, mladin, 'cak')
    
    #age 5 - 10: pb
    for i in range((2 + cak), (2 + cak + pbUp)):
        ped.izberi_poEBV_top_age_naive(i, 4, 'pb')
    
    
    #FEMALES
    #age 0
    #določi ženska teleta pod 12
    ped.izberi_poEBV_top_age_naive(0, telFn, 'telF')
    
    #age1
    #določi plemenske telice
    ped.izberi_poEBV_top_age_naive(1, ptn, 'pt')
    
    #age2
    for i in range(2, (1 + bmOdbira)):
        ped.izberi_poEBV_top_age_naive(i, ptn, 'k')
    
    #age3,4,5
    #odberi plemenske bm najprej
    for i in range((1 + bmOdbira), (1 + bmOdbira + bmUp)):
        ped.izberi_poEBV_top_age_naive(i, int(bmn / bmUp), 'pBM')
        ped.izberi_poEBV_top_age_naive(i, (ptn - int(bmn / bmUp)), 'k')
    
    #age 6
    #izberi odslužene bm
    ped.izberi_poEBV_top_age_naive((1 + bmOdbira + bmUp), int(bmn / bmUp), 'bm')
    
    
    #ostali so izločeni
    #določi spol ženskim živalim
    ped.set_sex_list(ped.row_cat('telF'), "F")
    ped.set_sex_list(ped.row_cat('pt'), "F")
    ped.set_sex_list(ped.row_cat('k'), "F")
    ped.set_sex_list(ped.row_cat('pBM'), "F")
    ped.set_sex_list(ped.row_cat('bm'), "F")
    
    
    #določi spol moškim živalim
    ped.set_sex_list(ped.row_cat('vhlevljeni'), "M")
    ped.set_sex_list(ped.row_cat('telM'), "M")
    ped.set_sex_list(ped.row_cat('bik12'), "M")
    ped.set_sex_list(ped.row_cat('mladi'), "M")
    ped.set_sex_list(ped.row_cat('cak'), "M")
    ped.set_sex_list(ped.row_cat('pb'), "M")
    ped.set_sex_list(ped.row_cat('pripust1'), "M")
    ped.set_sex_list(ped.row_cat('pripust2'), "M")
    
    #določi še izločene
    ped.set_cat_list(ped.row_cat(""), 'izl')
    
    ped.add_new_gen_naive(stNB, potomciNPn*2)
    
    ped.compute_age()
    #dodaj matere
    doloci_matere(ped)
    #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()
    #dodaj očete
    doloci_ocete(ped)
    #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()   
    
    categories.clear()
    categories = ped.save_cat()
    sex = ped.save_sex()
    active = ped.save_active()
    
    ped.write_ped("/home/jana/bin/AlphaSim1.05Linux/Pedigree_Gen1_StartWithTotal.txt")

