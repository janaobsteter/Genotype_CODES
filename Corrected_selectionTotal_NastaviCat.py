# -*- coding: utf-8 -*-
def selekcija_total(pedFile, **kwargs):
    print kwargs
    ped = pedigree(pedFile)
    
    #tukaj potem pridobi kategorije - če imaš samo eno burn-in, štartaš iz nule
    if max(ped.gen) == 1:
        ped.set_cat_gen(max(ped.gen), "nr")  # to je samo na prvem loopu
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 == 0], "F")
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 != 0], "M")
        ped.izberi_poEBV_top_catCurrent("F", int(kwargs.get('potomciNPn')), 'nr', 'potomciNP')
        ped.izberi_poEBV_top_catCurrent("M", int(kwargs.get('potomciNPn')), 'nr', 'potomciNP')
       
        #global categories #to moraš dat global samo v prvenm loopu, drugje dobiš return
        categories = ped.save_cat()
        #global sex
        sex = ped.save_sex()
        active = ped.save_active()
        ped = pedigree(pedFile) 
 
    elif max(ped.gens()) > 1:
        categories = ped.create_categoriesDict('Categories_gen' + str(max(ped.gens())) + 'DF.csv')  
        sex = ped.create_sexDict('Sex_gen' + str(max(ped.gens())) + 'DF.csv')  
        active = ped.create_activeDict('Active_gen' + str(max(ped.gens())) + 'DF.csv')  
       
    ped.set_sex_prevGen(sex)  # add sex information for individuals from prevGen
    ped.set_active_prevGen(active)  # add sex information for individuals from prevGen

    #remove category information from the ped itself
    for i in ped.gens():
        ped.set_cat_gen(i, "")

    #transfer culled (izlocene) category from prevGen
    ped.set_cat_old('izl', 'izl', categories)

    #compute age of the animals in the current selection year
    ped.compute_age()
    
    #################################################
    # FEMALES
    #################################################
    # age 0 - here you have newborn females (NB & potomkeNP) --> nekaj jih izloči, druge gredo naprej do ženskih telet
    ped.set_cat_sex_old("F", "potomciNP", "telF", categories) #potomke načrtnih parjenj gredo v telice
    izlF = int(kwargs.get('nrFn')) - int(kwargs.get('telFn'))  # number of culles NB females
    ped.izberi_poEBV_top("F", kwargs.get('telFn'), "nr", "telF", categories)  # izberi NB ženske, ki preživijo in postanejo telice
    ped.izloci_poEBV("F", izlF, "nr", categories)  # cull females (lowest on EBV) tukaj jih izloči, funkcija v modulu

    # age 1 - pri enem letu osemeni določeno število telic (% določen zgoraj), druge izloči
    if 'telF' in categories.keys():
        ped.izberi_poEBV_top("F", kwargs.get('ptn'), 'telF', 'pt', categories) #plemenske telice
        ped.izloci_poEBV("F", (len(categories['telF']) - kwargs.get('ptn')), 'telF', categories)  #preostale izloči

    # age > 2 - tukaj odbiraš in izločaš krave, odbiraš in izločaš BM
    # najprej dodaj nove krave, če jih že imaš v populaciji
    if ('pt' in categories.keys()): #če imaš v pedigreju plemenske telice
        ped.set_cat_old('pt', 'k', categories)  # osemenjene telice postanejo krave - predpostavimo, da vse
    # krave po 1., 2., 3. laktaciji prestavi naprej v krave - OZIROMA PODALJŠAJ STATUS!
    for i in range(2 + 1, (2 + kwargs.get('kraveUp'))):  # 2 + 1 - pri dveh letih prva laktacija, prestavljati začneš leto po tem
        ped.set_cat_age_old(i, 'k', 'k', categories)
    # potem izloči najstarejše krave - po 4. laktaciji
    if ('k' in categories.keys()) and ((kwargs.get('kraveUp') + 2) in ped.age()):  # izloči koliko laktacij + 2 leti
        ped.izloci_age_cat((kwargs.get('kraveUp') + 2), 'k', categories)


    # če imaš že dovolj stare krave, potem odberi BM
    # BM se odbira po drugi laktaciji - to je starost 3 - 4 (starost v pedigreju = 3, ker imaš tudi 0)
    if ('k' in categories.keys()) and ((1 + kwargs.get('bmOdbira')) in ped.age()):
        ped.izberi_poEBV_top_age("F", kwargs.get('bmOdbira')+1, int(kwargs.get('bmn') / kwargs.get('bmUp')), 'k', 'pBM', categories)  # izberi BM, ki jih osemeniš (plemenske BM = pBM) iz krav po 2. laktaciji
    # in izloči najastarejše BM, če jih imaš
    if ('bm' in categories.keys()):
        ped.izloci_cat('bm', categories)
    # ostale BM prestavi naprej - BM po 1. do izločitvene laktacije
    if 'pBM' in categories.keys():
        for i in range((1 + kwargs.get('bmOdbira') + 1), (
                1 + kwargs.get('bmOdbira') + kwargs.get('bmUp'))):  # 1 leto prva osemenitev, bm odbrane po 2. laktaciji, +1 da začneš prestavljat
            ped.set_cat_age_old(i, 'pBM', 'pBM', categories)
        # spremeni kategorijo iz plemenskih BM v bm v zadnji laktaciji 
        ped.set_cat_age_old((1 + kwargs.get('bmOdbira') + kwargs.get('bmUp')), 'pBM', 'bm',
                            categories)  

    #################################################################
    # MALES
    #################################################################
    # age 0: štartaš z NB in potomci NP --> odbereš vhlevljene iz potomcev NP in moška teleta in NB
    #NAJPREJ DELI, KI SO SKUPNI PROGENEMU IN GENOMSKEMU TESTIRANJU
    ped.izberi_random("M", kwargs.get('telMn'), "nr", "telM", categories) #izberi moška teleta, ki preživijo (random)
    ped.izloci_random("M", int(kwargs.get('nrMn') - kwargs.get('potomciNPn')-kwargs.get('telMn')), "nr", categories) #druga teleta izloči
    
    if 'telM' in categories.keys():
        ped.izberi_random("M", kwargs.get('bik12n'), 'telM', 'bik12', categories) #random odberi bike, ki preživijo do 2. leta
        ped.izloci_random("M", (len(categories['telM']) - kwargs.get('bik12n')), 'telM', categories) #izloči preostale
        
    if 'bik12' in categories.keys(): #izloci bike nad 2. leti
        ped.izloci_cat('bik12', categories)
    
    #Tukaj deli, kjer se progeni in genomski testi razlikujeta
    #progeni: potomciNP -> vhlevljeni -> mladi -> čakajoči -> pozitivno testirani
    #če je PROGENI TESTIRANJE
    if kwargs.get('EBV'):
        ped.izberi_poEBV_top("M", kwargs.get('vhlevljenin'), "potomciNP", "vhlevljeni",
                         categories)  # vhlevi najboljše potomceNP
        ped.izloci_poEBV("M", int(kwargs.get('potomciNPn') - kwargs.get('vhlevljenin')), 'potomciNP', categories) #druge potomceNP izloči

    
        # age1: tukaj odbereš mlade iz vhlevljenih bikov in bike, ki preživijo do drugega leta
        if 'vhlevljeni' in categories.keys(): #če imaš vhlevljene bike (samo v progenem testu)
            ped.izberi_poEBV_top("M", kwargs.get('mladin'), "vhlevljeni", "mladi", categories)  # odberi mlade
            ped.izberi_poEBV_OdDo("M", kwargs.get('mladin'), kwargs.get('vhlevljenin'), "vhlevljeni", "pripust1", categories)  # preostali vhlevljeni gredo v pripust

    
        # age > 2: tukaj mladi biki postanejo cakajoci in cakajo v testu
        #po koncanem testu odberes pozitivno testirane PB
        # mladi biki postanejo čakajoči (~1 leto, da se osemeni krave s semenom oz. malo po 2. letu)
        if 'mladi' in categories.keys():
            ped.set_cat_old('mladi', 'cak', categories) #mlade prestavi v cakajoce in jih izloci iz populacije
            ped.set_active_cat('mladi', 2, categories)

    
        # čakajočim bikov podaljšaj status (do starosti 5 let oz. kolikor let v testu)
        # hkrati jim tudi nastavi status izl
        # ped.set_cat_age_old(2, 'cak', 'cak', categories)
        if 'cak' in categories.keys():
            for i in range((2 + 1), (2 + kwargs.get('cak'))):  # 1 leto, ko začnejo semenit in so mladi biki, 3 so čakajoči, +1 da začneš prestavlajt
                ped.set_cat_age_old(i, 'cak', 'cak', categories)
    
    
    
        # če že imaš bike dovolj dolgo v testu, odberi pozitivno testirane bike
        if ('cak' in categories.keys()) and ((kwargs.get('cak') + 2) in ped.age()):  # +2 - eno leto so teleta, eno leto mladi biki
            ped.izberi_poEBV_top_age("M", (kwargs.get('cak') + 2), kwargs.get('pbn'), 'cak', 'pb', categories)
            ped.set_active_cat('cak', 2,
                            categories)  # tukaj moraš to nastaviti, zato ker fja izberi avtomatsko nastavi na active=1, vsi cakajoci so izloceni
            ped.izloci_poEBV_age("M", (kwargs.get('cak') + 2), kwargs.get('pbn'), 'cak', categories)  # TUKAJ MORA BITI ŠE STAROST!!!!!!!!!!!
    
    #genomski test: potomciNP = genomsko testiranje -> pozitivno testirani 
    if kwargs.get('gEBV'): #v prvem letu so vsi potomciNP v genomskem testiranju oz. pridobivanju gEBV
        ped.set_cat_sex_old('M', "potomciNP", "genTest", categories)     
        
        if 'genTest' in categories.keys(): #če imaš genomsko testirane bike 
            ped.izberi_poEBV_top("M", kwargs.get('pbn'), "genTest", "gpb", categories)  # odberi genomsko testirane bike za AI
            ped.izberi_poEBV_OdDo("M", kwargs.get('pbn'), kwargs.get('potomciNPn'), "genTest", "pripust1", categories)  # preostali vhlevljeni gredo v pripust

            

    #pripust in pb so spet enaki pri obeh testiranjih
    # povprečna doba v pripustu - glede na to odberi bike, ki preživijo še eno leto
    if 'pripust1' in categories.keys():
        ped.izberi_random("M", kwargs.get('pripust2n'), 'pripust1', 'pripust2', categories) #prestavi v 2. leto pripusta (ne vse - % glede na leta v UP)
        ped.izloci_random("M", (kwargs.get('pripust1n') - kwargs.get('pripust2n')), 'pripust1', categories) #preostale iz pripusta izloci

    if 'pripust2' in categories.keys(): #izloci po 2. letu v pripustu
        ped.izloci_cat('pripust2', categories)
        
    # plemenske bike prestavljaj naprej - zvseskozi, za očete pa nato uporabi le tiste iz željenih let
    if 'pb' in categories.keys():
        ped.set_cat_old('pb', 'pb', categories)
    if 'gpb' in categories.keys():
        ped.set_cat_old('gpb', 'gpb', categories)
    
    
    print ped.cat()
    #########################################################
    #add new generation
    #########################################################
    #tukaj potem dodaj eno generacijo novorojenih    
    ped.add_new_gen_naive(kwargs.get('stNBn'), kwargs.get('potomciNPn')*2)
    #določi starost glede na te novorojene
    ped.compute_age()
    #dodaj matere
    ped.doloci_matere(kwargs.get('stNBn'), kwargs.get('ptn'), kwargs.get('kraveUp'))
    #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()
    #dodaj očete
    ped.doloci_ocete(kwargs.get('stNBn'), kwargs.get('potomciNPn'), kwargs.get('cak'), kwargs.get('pripustDoz'), kwargs.get('pbUp'), kwargs.get('mladiDoz'), kwargs.get('pozitivnoTestDoz'))

    #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()   

    categories.clear()
    ped.save_cat_DF()
    ped.save_sex_DF()
    ped.save_active_DF()
    ped.write_pedTotal('/home/jana/PedTotal.txt')#("/home/jana/bin/AlphaSim1.05Linux/ExternalPedigree.txt")
    
    return ped, ped.save_cat(), ped.save_sex(), ped.save_active()

 #################

    
    
    
    
def nastavi_cat (PedFile, **kwargs):
    ped = pedigree(PedFile)
    ped.compute_age()
    
    
    #MALES FIRST
    #najprej deli, ki so skupni progenemu testi in genomskemu testiranju - to je neselekcionirana populacija
    #določi moška teleta pod 12
    ped.izberi_random_age_naive(0, kwargs.get('telMn'), 'telM')
    #določi bike nad 12 m
    ped.izberi_random_age_naive(1, kwargs.get('bik12n'), 'bik12')
    
    #PROGENI TEST
    if kwargs.get('EBV'):
        #age 0
        ped.izberi_poEBV_top_age_naive(0, kwargs.get('vhlevljenin'), 'vhlevljeni')     #določi vhlevljene
        #age1
        #določi mlade
        ped.izberi_poEBV_top_age_naive(1, kwargs.get('mladin'), 'mladi')
        #določi pripust - 1. leto
        ped.izberi_poEBV_OdDo_age_naive(1, kwargs.get('mladin'), kwargs.get('vhlevljenin'), 'pripust1')
        #age2,3,4
        for i in range(2, 2+kwargs.get('cak')): #leta, ko so cakajoci
            ped.izberi_poEBV_top_age_naive(i, kwargs.get('mladin'), 'cak')
        
        #od 1-2 leta v pripustu
        ped.izberi_poEBV_OdDo_age_naive(2, kwargs.get('mladin'), (kwargs.get('mladin') + kwargs.get('pripust2n')), 'pripust2')
        
        #age 5 - 10: pb
        pbAge = range((2 + kwargs.get('cak')), (2 + kwargs.get('cak')+ kwargs.get('pbUp'))) if (2 + kwargs.get('cak')+ kwargs.get('pbUp')) <= max(ped.gens()) else range((2 + kwargs.get('cak')), max(ped.gens()))
        for i in pbAge:
            ped.izberi_poEBV_top_age_naive(i, 4, 'pb')
    
    if kwargs.get('gEBV'):
        ped.izberi_poEBV_top_age_naive(0, kwargs.get('potomciNPn'), 'genTest')   
        
        #določi pripust - 1. leto
        ped.izberi_poEBV_OdDo_age_naive(1, kwargs.get('pbn'), kwargs.get('potomciNPn'), 'pripust1') #kateri niso odbrani po genomskih, gredo za pripust
        
        #od 1-2 leta v pripustu
        ped.izberi_poEBV_OdDo_age_naive(2, kwargs.get('pbn'), (kwargs.get('pbn') + kwargs.get('pripust2n')), 'pripust2')
        
        for i in range(1, 1+kwargs.get('pbUp')):
            ped.izberi_poEBV_top_age_naive(i, kwargs.get('pbn'), 'gpb') #odberi genomsko testirane bike za AI
    
    #FEMALES
    #age 0
    #določi ženska teleta pod 12
    ped.izberi_poEBV_top_age_naive(0, kwargs.get('telFn'), 'telF')
    
    #age1
    #določi plemenske telice
    ped.izberi_poEBV_top_age_naive(1, kwargs.get('ptn'), 'pt')
    
    #age2
    for i in range(2, (1 + kwargs.get('bmOdbira'))):
        ped.izberi_poEBV_top_age_naive(i, kwargs.get('ptn'), 'k')
    
    #age3,4,5
    #odberi plemenske bm najprej
    for i in range((1 + kwargs.get('bmOdbira')), (1 + kwargs.get('bmOdbira')+ kwargs.get('bmUp'))):
        ped.izberi_poEBV_top_age_naive(i, int(kwargs.get('bmn') / kwargs.get('bmUp')), 'pBM')
        ped.izberi_poEBV_top_age_naive(i, (kwargs.get('ptn') - int(kwargs.get('bmn') / kwargs.get('bmUp'))), 'k')
    
    #age 6
    #izberi odslužene bm
    ped.izberi_poEBV_top_age_naive((1 + kwargs.get('bmOdbira') + kwargs.get('bmUp')), int(kwargs.get('bmn') / kwargs.get('bmUp')), 'bm')
    
    
    #ostali so izločeni
    #določi spol ženskim živalim
    ped.set_sex_list(ped.row_cat('telF'), "F")
    ped.set_sex_list(ped.row_cat('pt'), "F")
    ped.set_sex_list(ped.row_cat('k'), "F")
    ped.set_sex_list(ped.row_cat('pBM'), "F")
    ped.set_sex_list(ped.row_cat('bm'), "F")
    
    
    #določi spol moškim živalim
    ped.set_sex_list(ped.row_cat('vhlevljeni'), "M")
    ped.set_sex_list(ped.row_cat('genTest'), "M")
    ped.set_sex_list(ped.row_cat('telM'), "M")
    ped.set_sex_list(ped.row_cat('bik12'), "M")
    ped.set_sex_list(ped.row_cat('mladi'), "M")
    ped.set_sex_list(ped.row_cat('cak'), "M")
    ped.set_sex_list(ped.row_cat('pb'), "M")
    ped.set_sex_list(ped.row_cat('gpb'), "M")
    ped.set_sex_list(ped.row_cat('pripust1'), "M")
    ped.set_sex_list(ped.row_cat('pripust2'), "M")
    
    #določi še izločene
    ped.set_sex_list(ped.row_cat(""), "I")
    ped.set_active_list(ped.row_cat(""), 2)
    ped.set_cat_list(ped.row_cat(""), 'izl')

    ped.add_new_gen_naive(kwargs.get('stNBn'), kwargs.get('potomciNPn')*2)
    
    ped.compute_age()
    #dodaj matere    
    ped.doloci_matere(kwargs.get('stNBn'), kwargs.get('ptn'), kwargs.get('kraveUp'))
    #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()
    #dodaj očete
    ped.doloci_ocete(kwargs.get('stNBn'), kwargs.get('potomciNPn'), kwargs.get('cak'), kwargs.get('pripustDoz'), kwargs.get('pbUp'), kwargs.get('mladiDoz'), kwargs.get('pozitivnoTestDoz'))
 
    ped.save_cat_DF()
    ped.save_sex_DF()
    ped.save_active_DF()
    ped.write_pedTotal('/home/jana/PedTotal.txt')#("/home/jana/bin/AlphaSim1.05Linux/ExternalPedigree.txt")
    
    return ped, ped.save_cat(), ped.save_sex(), ped.save_active()